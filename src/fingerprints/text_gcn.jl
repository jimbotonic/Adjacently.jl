#
# Text-GCN benchmark loader (R8 / R52 / 20NG / Ohsumed / MR).
#
# Reads the pre-tokenized corpus from Yao et al. 2019
# (github.com/yao8839836/text_gcn). Each dataset has two parallel files:
#
#   data/<DATASET>.txt          — one line per doc: "<id>\t<train|test>\t<label>"
#   data/corpus/<DATASET>.clean.txt — cleaned text, one doc per line, parallel to above
#
# The train/test split is fixed (matches the published benchmark numbers).
# Returns a vector of TextGcnDoc records and an Int label dictionary so the
# rest of the pipeline (build_vocab, compute_ppmi, build_domain_graph, NDF
# training) can run unchanged.

"""
    TextGcnDoc

A single Text-GCN benchmark document.

Fields:
- `tokens` — pre-tokenized word list (whitespace-split from the clean.txt line)
- `label`  — class string (e.g. "earn", "acq" for R8)
- `split`  — `:train` or `:test`
"""
struct TextGcnDoc
    tokens::Vector{String}
    label::String
    split::Symbol
end

"""
    read_text_gcn_corpus(dataset::AbstractString;
                        data_dir::AbstractString = joinpath(@__DIR__, "..", "..",
                            "datasets", "text_gcn", "text_gcn", "data"),
                        min_tokens::Int = 1)

Load a Text-GCN benchmark by name (`"R8"`, `"R52"`, `"20ng"`, `"ohsumed"`,
`"mr"`). `data_dir` points at the `data/` subdir of the cloned Yao 2019
repo. Documents with fewer than `min_tokens` tokens are dropped.

Returns `(docs, label2id)` where `docs::Vector{TextGcnDoc}` is in original
file order (so the train/test boundary is contiguous), and `label2id` maps
class strings to dense integer ids 1..n_classes (sorted by label string for
determinism).
"""
function read_text_gcn_corpus(dataset::AbstractString;
                              data_dir::AbstractString = joinpath(@__DIR__, "..", "..",
                                  "datasets", "text_gcn", "text_gcn", "data"),
                              min_tokens::Int = 1)
    split_path = joinpath(data_dir, "$dataset.txt")
    text_path  = joinpath(data_dir, "corpus", "$dataset.clean.txt")
    isfile(split_path) || error("missing split file: $split_path")
    isfile(text_path)  || error("missing corpus file: $text_path")

    splits = readlines(split_path)
    texts  = readlines(text_path)
    length(splits) == length(texts) ||
        error("length mismatch: $(length(splits)) split rows vs $(length(texts)) text rows")

    docs = TextGcnDoc[]
    labels = Set{String}()
    for (line, text) in zip(splits, texts)
        cols = split(line, '\t')
        length(cols) >= 3 || continue
        # cols = [doc_id, split, label]
        # Split tags vary across the Yao 2019 datasets:
        #   R8/R52/MR   → "train" / "test"
        #   Ohsumed     → "training" / "test"
        #   20NG        → "20news-bydate-train" / "20news-bydate-test"
        # Match permissively: anything containing "train" is train; "test" is test.
        tag = lowercase(cols[2])
        split_sym = occursin("train", tag) ? :train :
                    occursin("test",  tag) ? :test  :
                    error("unknown split tag: $(cols[2])")
        label = String(cols[3])
        tokens = split(text)
        length(tokens) >= min_tokens || continue
        push!(docs, TextGcnDoc(String.(tokens), label, split_sym))
        push!(labels, label)
    end
    # Sort labels for a deterministic id mapping.
    sorted_labels = sort!(collect(labels))
    label2id = Dict{String,Int}(l => i for (i, l) in enumerate(sorted_labels))
    return docs, label2id
end

"""
    textgcn_to_documents(docs, vocab, label2id; weight=:freq)

Convert `Vector{TextGcnDoc}` → `(train_docs, test_docs, n_classes)` where each
inner vector is `Vector{Document}` ready for NDF training. `weight` controls
the seed weighting in v_k:

- `:freq`    — normalized frequency vector v_k[u] = count(u in doc) / |doc|
- `:binary`  — v_k[u] = 1/|seed set| for words present (matches v1 DF semantics)
"""
function textgcn_to_documents(docs::Vector{TextGcnDoc},
                              vocab::Dict{String,Int},
                              label2id::Dict{String,Int};
                              weight::Symbol = :freq)
    weight ∈ (:freq, :binary) ||
        throw(ArgumentError("weight must be :freq or :binary"))
    train_docs = Document[]
    test_docs  = Document[]
    for d in docs
        # Map tokens to vertex ids, dropping OOV (only words in vocab survive).
        seed_counts = Dict{Int,Float32}()
        for w in d.tokens
            id = get(vocab, w, 0)
            id == 0 && continue
            seed_counts[id] = get(seed_counts, id, 0.0f0) + 1.0f0
        end
        isempty(seed_counts) && continue
        total = sum(values(seed_counts))
        seed_ids = collect(keys(seed_counts))
        seed_weights = if weight === :freq
            Float32[seed_counts[i] / total for i in seed_ids]
        else  # :binary
            Float32[1.0f0 / length(seed_ids) for _ in seed_ids]
        end
        doc = Document(seed_ids, seed_weights, label2id[d.label])
        if d.split === :train
            push!(train_docs, doc)
        else
            push!(test_docs, doc)
        end
    end
    return train_docs, test_docs, length(label2id)
end

"""
    load_bert_word_features(path, vocab) -> Matrix{Float32}

Load BERT averaged-contextual word embeddings from a `.npz` file produced by
`scripts/extract_bert_word_embeddings.py` and align them to `vocab`.

The script writes two parallel files:
- `<base>.npz`         containing key `X` : float32 (V_bert × H) embeddings
- `<base>.vocab.txt`   one token per line, parallel to rows of X

(We split vocab out of the .npz because NPZ.jl can't read NumPy's `object`
dtype string arrays.)

Returns a `(length(vocab), H)` matrix where row `vocab[w]` holds the BERT
vector for token `w`. Tokens missing from the BERT vocab (e.g. min_freq
mismatch) get a zero row; we warn so the caller can spot misalignment.
"""
function load_bert_word_features(path::AbstractString, vocab::Dict{String,Int})
    isfile(path) || error("missing BERT embedding file: $path")
    vocab_path = replace(path, r"\.npz$" => ".vocab.txt")
    isfile(vocab_path) ||
        error("missing parallel vocab file: $vocab_path (expected next to $path)")
    data = NPZ.npzread(path)
    haskey(data, "X") ||
        error("BERT npz must contain key 'X'; got $(collect(keys(data)))")
    X_bert = data["X"]
    bert_vocab = String[strip(line) for line in eachline(vocab_path) if !isempty(strip(line))]
    size(X_bert, 1) == length(bert_vocab) ||
        error("BERT X has $(size(X_bert,1)) rows but vocab has $(length(bert_vocab))")
    bert_to_idx = Dict{String,Int}(w => i for (i, w) in enumerate(bert_vocab))
    H = size(X_bert, 2)
    X = zeros(Float32, length(vocab), H)
    n_hit = 0
    for (w, vid) in vocab
        bidx = get(bert_to_idx, w, 0)
        bidx == 0 && continue
        X[vid, :] = X_bert[bidx, :]
        n_hit += 1
    end
    coverage = n_hit / length(vocab)
    coverage < 0.95 &&
        @warn "Low BERT vocab coverage" hit=n_hit vocab=length(vocab) coverage=round(coverage; digits=3) path=path
    return X
end

"""
    build_text_gcn_vocab(docs; min_freq=5)

Build a deterministic vocab dict (token → 1..V) from `docs`, keeping only
tokens that appear at least `min_freq` times **in the training split**.
Test-set tokens that don't appear in the training corpus are dropped (OOV)
to avoid label leakage. Sorting follows the same (`-count`, `token`) rule
as `build_vocab` so the vocab is reproducible across runs.
"""
function build_text_gcn_vocab(docs::Vector{TextGcnDoc}; min_freq::Int = 5)
    counts = Dict{String,Int}()
    for d in docs
        d.split === :train || continue
        for w in d.tokens
            counts[w] = get(counts, w, 0) + 1
        end
    end
    survivors = [(w, c) for (w, c) in counts if c >= min_freq]
    sort!(survivors; by = x -> (-x[2], x[1]))
    return Dict{String,Int}(w => i for (i, (w, _)) in enumerate(survivors))
end

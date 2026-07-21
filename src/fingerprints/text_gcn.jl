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
    to_char_ngrams(doc::TextGcnDoc; nmin=3, nmax=5) -> TextGcnDoc

Copy of `doc` with word tokens replaced by char-n-gram tokens (see
`char_ngram_tokens` / GLOSSARY.md).
"""
to_char_ngrams(d::TextGcnDoc; nmin::Int=3, nmax::Int=5) =
    TextGcnDoc(char_ngram_tokens(d.tokens; nmin=nmin, nmax=nmax), d.label, d.split)

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
- `:tfidf`   — TF-IDF weights normalized to sum to 1. Document frequency is
  computed from the **training split only** to avoid test-set leakage. Rare
  class-specific words get higher mass; generic words are down-weighted.
  Recommended for multi-class datasets with class-distinctive vocabulary
  (R52, 20NG, Ohsumed).
"""
function textgcn_to_documents(docs::Vector{TextGcnDoc},
                              vocab::Dict{String,Int},
                              label2id::Dict{String,Int};
                              weight::Symbol = :freq)
    weight ∈ (:freq, :binary, :tfidf, :tfidf_l2, :binary_raw) ||
        throw(ArgumentError("weight must be :freq, :binary, :tfidf, :tfidf_l2, or :binary_raw"))

    # Precompute IDF from the training split only (no test leakage).
    idf = if weight === :tfidf || weight === :tfidf_l2
        df = zeros(Int, length(vocab))
        n_train_docs = 0
        for d in docs
            d.split === :train || continue
            n_train_docs += 1
            seen = Set{Int}()
            for w in d.tokens
                id = get(vocab, w, 0)
                id == 0 && continue
                if !(id in seen)
                    push!(seen, id)
                    df[id] += 1
                end
            end
        end
        # Standard smoothed IDF: log((N + 1) / (df + 1)) + 1
        Float32[log(Float32(n_train_docs + 1) / (Float32(df[i]) + 1f0)) + 1f0
                for i in 1:length(vocab)]
    else
        Float32[]
    end

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
        elseif weight === :tfidf
            # TF-IDF, normalized to sum to 1 so seed vector remains a
            # probability distribution (preserves v_k semantics).
            raw = Float32[(seed_counts[i] / total) * idf[i] for i in seed_ids]
            s = sum(raw)
            s > 0f0 ? raw ./ s : Float32[1.0f0 / length(seed_ids) for _ in seed_ids]
        elseif weight === :tfidf_l2
            # L2-normalized TF-IDF — matches WideMLP's feature geometry (Galke &
            # Scherp 2022) for a fair-input comparison. Breaks the sum-to-1 v_k
            # probability semantics on purpose; keep :tfidf/:freq for v1-faithful.
            raw = Float32[(seed_counts[i] / total) * idf[i] for i in seed_ids]
            nrm = sqrt(sum(abs2, raw))
            nrm > 0f0 ? raw ./ nrm : Float32[1.0f0 / length(seed_ids) for _ in seed_ids]
        elseif weight === :binary_raw
            # Raw presence (unnormalized 1s) — WideMLP-style binary BoW.
            ones(Float32, length(seed_ids))
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
    read_blog_posts(dir; min_post_tokens=10, n_authors=500,
                    min_posts_per_author=20, test_frac=0.2,
                    chunk_words=0, label_by=:author,
                    stem=false, seed=0) -> Vector{TextGcnDoc}

Layout for per-post or per-chunk classification on the blog corpus.
Returns one TextGcnDoc per "document" labelled by `blogger_id`
(authorship attribution, the v1 §4.2 task) or by `gender` (binary
gender classification on chunks). Pipeline:

1. Walk every XML, extract each `<post>` as its own raw chunk; drop posts
   with fewer than `min_post_tokens` tokens.
2. Filter to authors with at least `min_posts_per_author` qualifying posts.
3. Select the top `n_authors` by post count (most data per class first).
4. **Document granularity** (controlled by `chunk_words`):
   - `chunk_words == 0` (default) — one document per individual post.
     Matches the v1 §4.2 setup but with per-post training.
   - `chunk_words > 0` — concatenate each author's posts in order, then
     split into successive ~`chunk_words`-token chunks. Each chunk is one
     document. Closer to the Koppel et al. (2011) setup where each test
     instance is a 500-word "snippet" rather than a single post. The last
     chunk may be shorter; dropped if shorter than `chunk_words / 4`.
5. For each selected author: shuffle documents with `MersenneTwister(seed)`,
   split `test_frac` to `:test`, rest to `:train`.

Returned vector is in (per-author, then per-doc) order. Caller typically
passes it to `build_text_gcn_vocab` and `textgcn_to_documents`.
"""
function read_blog_posts(dir::AbstractString;
                         min_post_tokens::Int=10,
                         n_authors::Int=500,
                         min_posts_per_author::Int=20,
                         test_frac::Float64=0.2,
                         chunk_words::Int=0,
                         label_by::Symbol=:author,
                         stem::Bool=false,
                         seed::Int=0)
    label_by ∈ (:author, :gender) ||
        throw(ArgumentError("label_by must be :author or :gender, got $label_by"))
    files = filter(f -> endswith(f, ".xml"), readdir(dir; join=true, sort=true))
    # author_id → Vector{post_tokens}, plus gender lookup
    by_author = Dict{String,Vector{Vector{String}}}()
    author_gender = Dict{String,Symbol}()
    for f in files
        try
            aid, gender, posts = _parse_blog_posts(f; stem=stem)
            keep = [p for p in posts if length(p) >= min_post_tokens]
            isempty(keep) && continue
            append!(get!(by_author, aid, Vector{Vector{String}}()), keep)
            author_gender[aid] = gender
        catch err
            @warn "Skipping $f: $err"
        end
    end

    # Filter authors with enough posts, take top N by count.
    eligible = [(aid, posts) for (aid, posts) in by_author
                if length(posts) >= min_posts_per_author]
    sort!(eligible; by = x -> (-length(x[2]), x[1]))
    n_keep = min(n_authors, length(eligible))
    selected = eligible[1:n_keep]

    # Per-author train/test split with a shared RNG (re-runs stable given
    # same seed; selection order above is already deterministic).
    rng = MersenneTwister(seed)
    out = TextGcnDoc[]
    for (aid, posts) in selected
        # Build the per-author document list, either one-doc-per-post (default)
        # or coalesce posts into ~chunk_words-token chunks.
        docs = if chunk_words <= 0
            posts
        else
            all_toks = reduce(vcat, posts)
            n = length(all_toks)
            chunks = Vector{Vector{String}}()
            i = 1
            while i <= n
                stop = min(n, i + chunk_words - 1)
                push!(chunks, all_toks[i:stop])
                i = stop + 1
            end
            # Drop the trailing chunk if it's much shorter than the target
            # — keeps per-doc size roughly uniform.
            if !isempty(chunks) && length(chunks[end]) < chunk_words ÷ 4 && length(chunks) > 1
                pop!(chunks)
            end
            chunks
        end
        isempty(docs) && continue
        perm = Random.shuffle(rng, collect(1:length(docs)))
        n_test = max(1, round(Int, test_frac * length(docs)))
        test_idx = Set(perm[1:n_test])
        label = label_by === :author ? aid : String(author_gender[aid])
        for (i, tokens) in enumerate(docs)
            split_sym = i in test_idx ? :test : :train
            push!(out, TextGcnDoc(tokens, label, split_sym))
        end
    end
    return out
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

"""
    word_class_weights(docs, vocab, label2id; smooth=1f0, floor=0.05f0) -> Vector{Float32}

Per-word class-discriminativeness, computed from the TRAINING split only (no
leakage). For each word `w`, take its class distribution `p(c|w)` over training
docs that contain `w`, and score it by the KL divergence from the class prior
`p(c)` — words concentrated in a few classes (e.g. class-specific jargon) score
high, stopword-like words spread across classes score ~0. Returned weights are
rescaled to `[floor, 1]` so no edge is fully zeroed. Used to reweight PPMI edges
so diffusion flows along class-consistent co-occurrences (class-aware graph).
"""
function word_class_weights(docs::Vector{TextGcnDoc}, vocab::Dict{String,Int},
                            label2id::Dict{String,Int};
                            smooth::Float32=1f0, floor::Float32=0.05f0)
    V = length(vocab); C = length(label2id)
    wc = zeros(Float32, V, C)          # wc[w,c] = # train docs of class c with word w
    classtot = zeros(Float32, C)
    for d in docs
        d.split === :train || continue
        c = label2id[d.label]
        classtot[c] += 1f0
        seen = Set{Int}()
        for t in d.tokens
            id = get(vocab, t, 0); id == 0 && continue
            if !(id in seen); push!(seen, id); wc[id, c] += 1f0; end
        end
    end
    ntrain = sum(classtot)
    pc = (classtot .+ smooth) ./ (ntrain + C * smooth)   # class prior
    score = Vector{Float32}(undef, V)
    @inbounds for i in 1:V
        row = @view wc[i, :]
        s = sum(row) + C * smooth
        kl = 0f0
        for c in 1:C
            pcw = (row[c] + smooth) / s
            kl += pcw * log(pcw / pc[c])
        end
        score[i] = kl
    end
    mx = maximum(score); mx = mx > 0f0 ? mx : 1f0
    return floor .+ (1f0 - floor) .* (score ./ mx)
end

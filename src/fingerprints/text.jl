#
# Adjacently: Julia Complex Directed Networks Library
# Copyright (C) 2016-2026 Jimmy Dubuisson <jimmy@dubuisson.ch>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#

# Text-pipeline helpers for the v1 DF blog experiments and the Text-GCN
# benchmark suite. Kept inside the Fingerprints module so the model and the
# pipeline that feeds it ship together.

"""
    BlogDoc(blogger_id, gender, age, industry, sign, tokens)

A single blogger's complete corpus (all posts concatenated) with metadata
parsed from the Schler-corpus filename and the post bodies tokenized.
"""
struct BlogDoc
    blogger_id::String
    gender::Symbol      # :male, :female
    age::Int
    industry::String
    sign::String
    tokens::Vector{String}
end

# Compact English stopword list (similar to NLTK's English stopwords minus the
# punctuation tokens we already strip). Kept inline to avoid a CSV dep.
const STOPWORDS = Set([
    "a","about","above","after","again","against","all","am","an","and","any",
    "are","as","at","be","because","been","before","being","below","between",
    "both","but","by","can","cannot","could","did","do","does","doing","don",
    "down","during","each","few","for","from","further","had","has","have",
    "having","he","her","here","hers","herself","him","himself","his","how",
    "i","if","in","into","is","it","its","itself","just","ll","m","me","more",
    "most","my","myself","no","nor","not","now","o","of","off","on","once",
    "only","or","other","our","ours","ourselves","out","over","own","re","s",
    "same","she","should","so","some","such","t","than","that","the","their",
    "theirs","them","themselves","then","there","these","they","this","those",
    "through","to","too","under","until","up","ve","very","was","we","were",
    "what","when","where","which","while","who","whom","why","will","with",
    "would","you","your","yours","yourself","yourselves",
])

"""
    simple_tokenize(text; remove_stopwords=true, min_token_length=2, stem=false) -> Vector{String}

Lowercase, split on non-alphanumeric, drop tokens shorter than
`min_token_length` characters, drop stopwords if requested. When `stem=true`,
apply a minimal Porter-style stemmer (`porter_stem_lite`) — sufficient to
collapse plural/tense/derivational variants for the v1 reproduction without
requiring an external Snowball dependency.
"""
function simple_tokenize(text::AbstractString;
                         remove_stopwords::Bool=true,
                         min_token_length::Int=2,
                         stem::Bool=false)
    tokens = String[]
    for tok in eachmatch(r"[A-Za-z]+", lowercase(text))
        s = tok.match
        length(s) < min_token_length && continue
        remove_stopwords && s in STOPWORDS && continue
        push!(tokens, stem ? porter_stem_lite(s) : s)
    end
    return tokens
end

# Minimal Porter-style stemmer. Covers the high-frequency English suffixes
# (plurals, tense, common derivational endings) without claiming full Porter
# compliance. Enough to give the v1 reproduction a similar vocab size to
# Schler's setup; falls short of full Porter on rare edge cases.
const _STEM_SUFFIXES_LONG = String[
    "ational", "tional", "fulness", "ousness", "iveness",
    "ization", "ational", "isation", "ization",
    "ousness", "fulness", "iveness", "biliti",
    "ation",  "ition",   "tions",   "ators",
    "ables",  "ation",   "ition",   "ness",
    "less",   "ment",    "able",    "ible",
    "ours",   "ators",   "ation",
]
const _STEM_SUFFIXES_SHORT = String[
    "ing", "ies", "ied", "edly", "ed", "ly", "es", "s",
]

function _strip_suffix(s::AbstractString, suffix::AbstractString)
    return s[1:end - length(suffix)]
end

"""
    porter_stem_lite(word) -> String

Reduce a word to a stem by stripping common English suffixes. Not full
Porter (no consonant-vowel measure, no double-letter rules) — a deliberate
~50-line approximation that covers the suffixes that matter most for
vocabulary collapse on the blog corpus.
"""
function porter_stem_lite(word::AbstractString)
    s = String(word)
    length(s) <= 3 && return s

    # Try long suffixes first (each consumes more characters).
    for suf in _STEM_SUFFIXES_LONG
        if length(s) > length(suf) + 2 && endswith(s, suf)
            s = _strip_suffix(s, suf)
            break
        end
    end

    # Then short suffixes; iterate until none applies (handles e.g. "studies" → "studi" → "studi").
    changed = true
    while changed && length(s) > 3
        changed = false
        for suf in _STEM_SUFFIXES_SHORT
            if length(s) > length(suf) + 2 && endswith(s, suf)
                # "ies" / "ied" → "i" (so "studies" → "studi", "tried" → "tri").
                if suf == "ies" || suf == "ied"
                    s = _strip_suffix(s, suf) * "i"
                elseif suf == "ing" || suf == "ed"
                    s = _strip_suffix(s, suf)
                    # Collapse doubled consonants (e.g. "running" → "runn" → "run").
                    if length(s) >= 2 && s[end] == s[end-1] && s[end] ∉ ('s', 'l', 'z')
                        s = s[1:end-1]
                    end
                else
                    s = _strip_suffix(s, suf)
                end
                changed = true
                break
            end
        end
    end
    return s
end

# Parse one blog XML file. Filename format:
#   <id>.<gender>.<age>.<industry>.<sign>.xml
# Body format: alternating <date>...</date> and <post>...</post> nodes wrapped
# in <Blog>...</Blog>. We extract just the post bodies via regex — the corpus
# XML is non-strict (entity errors, mixed encodings) so a strict parser
# breaks; regex is more forgiving.
function _parse_blog_file(path::AbstractString; stem::Bool=false)
    fname = basename(path)
    parts = split(replace(fname, r"\.xml$" => ""), ".")
    length(parts) >= 5 ||
        throw(ArgumentError("Unexpected blog filename: $fname"))

    blogger_id = parts[1]
    gender_str = lowercase(parts[2])
    gender = gender_str == "male" ? :male :
             gender_str == "female" ? :female :
             throw(ArgumentError("Unrecognized gender '$gender_str' in $fname"))
    age = parse(Int, parts[3])
    industry = parts[4]
    sign = parts[5]

    raw = read(path, String)
    # Strip BOM, normalize line endings.
    raw = replace(raw, "﻿" => "", "\r\n" => "\n", "\r" => "\n")

    tokens = String[]
    for m in eachmatch(r"<post>(.*?)</post>"s, raw)
        post_body = m.captures[1]
        append!(tokens, simple_tokenize(post_body; stem=stem))
    end

    return BlogDoc(blogger_id, gender, age, industry, sign, tokens)
end

# Parse a blog XML into individual posts (one tokenized string per <post>
# block). Used by the authorship-attribution loader which treats each post
# as a separate document with the blogger_id as label, in contrast to
# `_parse_blog_file` which concatenates all posts of a blogger.
# Returns (blogger_id, gender_symbol, posts) so callers can label by either.
function _parse_blog_posts(path::AbstractString; stem::Bool=false)
    fname = basename(path)
    parts = split(replace(fname, r"\.xml$" => ""), ".")
    length(parts) >= 5 ||
        throw(ArgumentError("Unexpected blog filename: $fname"))
    blogger_id = parts[1]
    gender_str = lowercase(parts[2])
    gender = gender_str == "male" ? :male :
             gender_str == "female" ? :female :
             throw(ArgumentError("Unrecognized gender '$gender_str' in $fname"))

    raw = read(path, String)
    raw = replace(raw, "﻿" => "", "\r\n" => "\n", "\r" => "\n")

    posts = Vector{Vector{String}}()
    for m in eachmatch(r"<post>(.*?)</post>"s, raw)
        toks = simple_tokenize(m.captures[1]; stem=stem)
        isempty(toks) && continue
        push!(posts, toks)
    end
    return blogger_id, gender, posts
end

"""
    read_blog_corpus(dir; max_docs=typemax(Int), min_tokens=8,
                          balance_genders=false, seed=0) -> Vector{BlogDoc}

Load blog files from `dir`. One file per blogger in Schler-corpus naming
format. Bloggers with fewer than `min_tokens` total tokens after tokenization
are dropped. If `balance_genders=true`, the returned vector contains an equal
number of male and female bloggers (random subsample with seed).
"""
function read_blog_corpus(dir::AbstractString;
                          max_docs::Int=typemax(Int),
                          min_tokens::Int=8,
                          balance_genders::Bool=false,
                          stem::Bool=false,
                          seed::Int=0)
    files = filter(f -> endswith(f, ".xml"), readdir(dir; join=true, sort=true))
    docs = BlogDoc[]
    for f in files
        try
            d = _parse_blog_file(f; stem=stem)
            length(d.tokens) < min_tokens && continue
            push!(docs, d)
        catch err
            @warn "Skipping $f: $err"
        end
        length(docs) >= max_docs && break
    end

    if balance_genders
        rng = MersenneTwister(seed)
        males = filter(d -> d.gender == :male, docs)
        females = filter(d -> d.gender == :female, docs)
        n = min(length(males), length(females))
        Random.shuffle!(rng, males); Random.shuffle!(rng, females)
        docs = vcat(males[1:n], females[1:n])
        Random.shuffle!(rng, docs)
    end
    return docs
end

"""
    build_vocab(docs; min_freq=5, max_vocab=typemax(Int)) -> Dict{String,Int}

Map tokens to 1-based vertex indices for the domain graph. Tokens appearing
fewer than `min_freq` times across the corpus are dropped, then the top
`max_vocab` by frequency are retained. The dict maps surviving tokens to
contiguous indices `1:|V|`.
"""
function build_vocab(docs::AbstractVector{BlogDoc};
                     min_freq::Int=5,
                     max_vocab::Int=typemax(Int))
    counts = Dict{String,Int}()
    for d in docs, t in d.tokens
        counts[t] = get(counts, t, 0) + 1
    end
    survivors = [(tok, c) for (tok, c) in counts if c >= min_freq]
    # Deterministic order: descending count, ties broken alphabetically by token.
    # Without the secondary key, Dict iteration order leaks into vertex IDs and
    # breaks checkpoint reuse across runs.
    sort!(survivors; by = x -> (-x[2], x[1]))
    length(survivors) > max_vocab && (survivors = survivors[1:max_vocab])
    return Dict(tok => i for (i, (tok, _)) in enumerate(survivors))
end

"""
    compute_ppmi(docs, vocab; window=20) -> SparseMatrixCSC{Float32}

Compute positive pointwise mutual information (PPMI) over co-occurrences
within a sliding window of `window` tokens. PPMI(u,v) = max(0, log(p(u,v) /
(p(u) p(v)))) using empirical windowed counts; this is the standard
Text-GCN-style edge weighting (Yao et al. 2019). Returns a symmetric sparse
matrix of size `|V| × |V|`.
"""
function compute_ppmi(docs::AbstractVector,
                      vocab::Dict{String,Int};
                      window::Int=20)
    V = length(vocab)
    pair_counts = Dict{Tuple{Int,Int},Int}()
    token_counts = zeros(Int, V)
    n_windows = 0

    for d in docs
        ids = Int[]
        for t in d.tokens
            id = get(vocab, t, 0)
            id != 0 && push!(ids, id)
        end
        n = length(ids)
        n == 0 && continue
        for i in 1:n
            stop = min(n, i + window - 1)
            in_window = unique(ids[i:stop])
            for u in in_window
                token_counts[u] += 1
            end
            for a in 1:length(in_window), b in a+1:length(in_window)
                u, v = in_window[a], in_window[b]
                u, v = min(u, v), max(u, v)
                pair_counts[(u, v)] = get(pair_counts, (u, v), 0) + 1
            end
            n_windows += 1
        end
    end

    # Build PPMI sparse matrix.
    rows = Int[]; cols = Int[]; vals = Float32[]
    for ((u, v), c) in pair_counts
        p_uv = c / n_windows
        p_u  = token_counts[u] / n_windows
        p_v  = token_counts[v] / n_windows
        pmi = log(p_uv / (p_u * p_v + 1.0f-12))
        pmi <= 0 && continue
        push!(rows, u); push!(cols, v); push!(vals, Float32(pmi))
        push!(rows, v); push!(cols, u); push!(vals, Float32(pmi))
    end
    return sparse(rows, cols, vals, V, V)
end

"""
    build_domain_graph(K::SparseMatrixCSC, density::Real; weighted::Bool=false)
        -> SparseMatrixCSC

Threshold a weighted association matrix `K` at the top-`density` fraction of
edges (per the v1 DF construction). With `weighted=false` (default, matches
the legacy behaviour) returns a 0/1 sparse matrix. With `weighted=true`
returns a sparse matrix preserving the original PPMI / collocation
magnitudes — this is faithful to the v1 paper, which explicitly allows a
weighted domain graph, and tends to help when the surviving edges have a
wide range of strengths (e.g. PPMI on multi-class datasets where rare-but-
strong co-occurrences dominate class structure).

Memory is the bottleneck for large `K` — the implementation avoids
`findnz`'s triple-copy and walks `K`'s internal CSC arrays directly,
allocating only the final `(keep_rows, keep_cols)` pair plus one transient
sort buffer. Ties at the threshold are kept (slight overshoot of the
target).
"""
function build_domain_graph(K::SparseMatrixCSC, density::Real;
                            weighted::Bool=false)
    n = size(K, 1)
    target = max(1, round(Int, density * n * (n - 1)))
    nzv = K.nzval
    n_entries = length(nzv)
    if n_entries <= target
        # Reuse K's structure; choose values per `weighted` flag.
        vals = weighted ? copy(nzv) : ones(Float32, n_entries)
        return SparseMatrixCSC(n, n, copy(K.colptr), copy(K.rowval), vals)
    end

    # Threshold = the target-th largest nzval. Copy + partial-sort the values
    # in-place on the copy, then release the copy before allocating keep arrays.
    sorted_vals = copy(nzv)
    partialsort!(sorted_vals, target; rev=true)
    threshold = sorted_vals[target]
    sorted_vals = nothing
    GC.gc()

    # Two-pass build to size the keep arrays exactly. Ties are included; on
    # heavy ties the actual edge count can exceed `target` by a small margin.
    keep_count = count(>=(threshold), nzv)
    keep_rows = Vector{Int}(undef, keep_count)
    keep_cols = Vector{Int}(undef, keep_count)
    keep_vals = weighted ? Vector{Float32}(undef, keep_count) : nothing
    pos = 1
    @inbounds for j in 1:n
        for idx in K.colptr[j]:(K.colptr[j + 1] - 1)
            v = nzv[idx]
            if v >= threshold
                keep_rows[pos] = K.rowval[idx]
                keep_cols[pos] = j
                if weighted
                    keep_vals[pos] = Float32(v)
                end
                pos += 1
            end
        end
    end
    final_vals = weighted ? keep_vals : ones(Float32, keep_count)
    return sparse(keep_rows, keep_cols, final_vals, n, n)
end

"""
    seed_from_doc(doc::BlogDoc, vocab) -> (indices, weights)

Build the sparse personalization vector v_k for a single document: vertex
indices and their normalized frequencies (sum to 1).
"""
function seed_from_doc(doc::BlogDoc, vocab::Dict{String,Int})
    counts = Dict{Int,Int}()
    total = 0
    for t in doc.tokens
        id = get(vocab, t, 0)
        id == 0 && continue
        counts[id] = get(counts, id, 0) + 1
        total += 1
    end
    total == 0 && return (Int[], Float32[])
    idx = sort!(collect(keys(counts)))
    w = Float32[counts[i] / total for i in idx]
    return idx, w
end

"""
    compute_v1_collocation(docs, vocab; β=1.0f0, σ=1.0f0, max_distance=10) -> SparseMatrixCSC{Float32}

Reproduce the v1 DF directed association matrix (Dubuisson et al. 2014, §2.1):

    K(k)_uv = g(h(u,v)) · ∑_{(i,j) ∈ s_uv}  f(i, j)

with f(i,j) = exp(-(j-i-1)^β / σ), g(x) = -log(x), and s_uv = positions of v
occurring strictly after each occurrence of u and strictly before u's next
occurrence in document k. h(u,v) = |s_uv| / Σ_{u',v'} |s_{u',v'}| is the
within-document frequency of the (u,v) pair. K(Σ) = ∑_k K(k) is returned as a
sparse, **directed** matrix.

`max_distance` caps the f-decay window (default 10 — at β=σ=1 the weight
falls below 1e-4 well before that). Reduces inner-loop cost from O(|T'|^2)
to O(|T'| · max_distance).
"""
function compute_v1_collocation(docs::AbstractVector{BlogDoc},
                                vocab::Dict{String,Int};
                                β::Float32=1.0f0,
                                σ::Float32=1.0f0,
                                max_distance::Int=10)
    V = length(vocab)
    K_total = Dict{Tuple{Int,Int},Float32}()

    for d in docs
        ids = Int[]
        for t in d.tokens
            id = get(vocab, t, 0)
            id != 0 && push!(ids, id)
        end
        n = length(ids)
        n < 2 && continue

        # Per-token occurrence positions.
        positions = Dict{Int,Vector{Int}}()
        for (i, id) in enumerate(ids)
            push!(get!(() -> Int[], positions, id), i)
        end

        # Per-doc accumulators.
        K_doc = Dict{Tuple{Int,Int},Float32}()
        s_doc = Dict{Tuple{Int,Int},Int}()

        for (u, u_pos) in positions
            # Sentinel n+1 marks "past end of doc".
            extended = vcat(u_pos, n + 1)
            for k in 1:length(u_pos)
                p_u_k     = u_pos[k]
                p_u_next  = extended[k + 1]
                # v's positions in (p_u_k, p_u_next), capped by max_distance.
                upper = min(n, p_u_k + max_distance, p_u_next - 1)
                for j in (p_u_k + 1):upper
                    v = ids[j]
                    v == u && continue
                    distance = j - p_u_k - 1
                    f_val    = exp(-(Float32(distance)^β) / σ)
                    key      = (u, v)
                    K_doc[key] = get(K_doc, key, 0.0f0) + f_val
                    s_doc[key] = get(s_doc, key, 0) + 1
                end
            end
        end

        total_s = sum(values(s_doc))
        total_s == 0 && continue

        # K(k)_uv = -log(h(u,v)) · K_doc[u,v], h = |s_uv| / total_s.
        for ((u, v), Kv) in K_doc
            h   = s_doc[(u, v)] / total_s
            g_h = -log(Float32(h))
            K_total[(u, v)] = get(K_total, (u, v), 0.0f0) + g_h * Kv
        end
    end

    rows = Int[]; cols = Int[]; vals = Float32[]
    for ((u, v), val) in K_total
        push!(rows, u); push!(cols, v); push!(vals, val)
    end
    return sparse(rows, cols, vals, V, V)
end

"""
    blogdocs_to_documents(docs, vocab; label_fn) -> Vector{Document}

Convert a list of `BlogDoc`s into `Document`s suitable for `train!`. The
`label_fn` callable maps `BlogDoc -> Int` and determines the class (1-based).
For gender detection: `label_fn = d -> d.gender == :male ? 1 : 2`.
"""
function blogdocs_to_documents(docs::AbstractVector{BlogDoc},
                               vocab::Dict{String,Int};
                               label_fn::Function)
    out = Document[]
    for d in docs
        idx, w = seed_from_doc(d, vocab)
        isempty(idx) && continue
        push!(out, Document(idx, w, label_fn(d)))
    end
    return out
end

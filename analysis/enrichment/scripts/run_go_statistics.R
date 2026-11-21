#!/usr/bin/env Rscript
# Frozen-protocol GO-BP statistics for the PRJNA734471 DSE manuscript.
#
# This script does not rerun RNA-seq models.  It accepts already defined
# sequence-equivalence units and direct, frozen host GO memberships only.
# Source-protein GO ledgers, candidate-only annotations and newly transferred
# annotations must NOT be used as input to this script.
#
# Primary questions (pre-specified):
#   A. Are the 391 opposite-direction units over-represented for direct BP
#      labels relative to the annotated non-risk universe?
#   B. Within the 391 units, do the two opposite-direction classes differ in
#      BP-label composition?
# The strict-127 analyses are nested sensitivity outputs only, never an
# independent validation or a second discovery screen.
#
# Usage (after protocol approval):
#   Rscript run_go_statistics.R \
#     --input-dir input --output-dir output \
#     --direction-a P0_pos_P1_neg --direction-b P0_neg_P1_pos
#
# Required input files in --input-dir:
#   universe.tsv    unit_id, baseMean, length, opposite, direction, strict127
#   memberships.tsv unit_id, term_id
#   terms.tsv       term_id, name, namespace
#   term_exclusions.tsv  term_id, reason; frozen non-biological BP labels to omit
#
# Statistical protocol locked in this file:
#   * direct GO-BP memberships only; no GO ancestor propagation;
#   * actual BP-annotated units define U, the testable background;
#   * roots excluded; term size is >=5 and <=95% of U;
#   * A: one-sided Fisher over-representation; B: two-sided Fisher;
#   * primary editorial decision is q_combined, BH over all eligible A and B
#     term tests; q_within_analysis is retained for transparency;
#   * strict-127 outputs contain all terms but make no new discovery calls;
#   * 10,000 term-wise, conditional-hypergeometric Monte Carlo draws
#     (seed 20261110), strata=length quintile x baseMean quintile.  A preserves
#     selected counts by strata; B exchanges direction labels inside selected
#     391 units while preserving stratum counts.  It is not a joint permutation
#     of the complete term-label matrix.

options(stringsAsFactors = FALSE, warn = 1)

parse_args <- function(x) {
  out <- list(
    input_dir = "input", output_dir = "output",
    direction_a = "P0up_P1down", direction_b = "P0down_P1up",
    n_perm = 10000L, seed = 20261110L, expected_universe = 41139L
  )
  i <- 1L
  while (i <= length(x)) {
    key <- x[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    key <- sub("^--", "", key)
    if (key %in% c("help", "h")) {
      cat("See header of run_go_statistics.R for usage.\n")
      quit(status = 0)
    }
    if (i == length(x)) stop("Missing value after --", key)
    value <- x[[i + 1L]]
    key <- gsub("-", "_", key)
    if (!(key %in% names(out))) stop("Unknown option --", gsub("_", "-", key))
    out[[key]] <- value
    i <- i + 2L
  }
  out$n_perm <- as.integer(out$n_perm)
  out$seed <- as.integer(out$seed)
  out$expected_universe <- as.integer(out$expected_universe)
  if (is.na(out$n_perm) || out$n_perm < 1000L) stop("--n-perm must be >=1000")
  if (is.na(out$seed)) stop("--seed must be an integer")
  if (is.na(out$expected_universe) || out$expected_universe < 1L) stop("--expected-universe must be a positive integer")
  out
}

cfg <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(cfg$output_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_required <- function(path, fields) {
  if (!file.exists(path)) stop("Required input missing: ", path)
  x <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
                  quote = "", comment.char = "", stringsAsFactors = FALSE)
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(basename(path), " lacks required column(s): ", paste(absent, collapse = ", "))
  x
}

as_flag <- function(x, field) {
  y <- trimws(tolower(as.character(x)))
  ok <- y %in% c("true", "t", "1", "yes", "y")
  no <- y %in% c("false", "f", "0", "no", "n", "") | is.na(y)
  if (any(!ok & !no)) stop("Unrecognised boolean value(s) in ", field, ": ",
                            paste(unique(y[!ok & !no]), collapse = ", "))
  ok
}

root_terms <- c("GO:0008150", "GO:0003674", "GO:0005575")

universe <- read_tsv_required(file.path(cfg$input_dir, "universe.tsv"),
                              c("unit_id", "baseMean", "length", "opposite", "direction", "strict127"))
memberships <- read_tsv_required(file.path(cfg$input_dir, "memberships.tsv"), c("unit_id", "term_id"))
terms <- read_tsv_required(file.path(cfg$input_dir, "terms.tsv"), c("term_id", "name", "namespace"))
term_exclusions <- read_tsv_required(file.path(cfg$input_dir, "term_exclusions.tsv"), c("term_id", "reason"))

if (anyDuplicated(universe$unit_id)) stop("universe.tsv has duplicated unit_id values")
if (anyDuplicated(terms$term_id)) stop("terms.tsv has duplicated term_id values")
if (any(is.na(universe$unit_id) | universe$unit_id == "")) stop("universe.tsv has blank unit_id")
if (nrow(universe) != cfg$expected_universe) {
  stop("universe.tsv has ", nrow(universe), " rows; frozen protocol expects ", cfg$expected_universe,
       ". Do not substitute a different universe without a revised protocol.")
}

universe$opposite <- as_flag(universe$opposite, "opposite")
universe$strict127 <- as_flag(universe$strict127, "strict127")
universe$baseMean <- suppressWarnings(as.numeric(universe$baseMean))
universe$length <- suppressWarnings(as.numeric(universe$length))

memberships <- unique(memberships[, c("unit_id", "term_id")])
memberships <- memberships[memberships$unit_id %in% universe$unit_id, , drop = FALSE]
memberships <- memberships[memberships$term_id %in% terms$term_id, , drop = FALSE]
terms$namespace_clean <- tolower(trimws(terms$namespace))
bp_terms <- terms[terms$namespace_clean %in% c("biological_process", "bp", "biological process") &
                    !(terms$term_id %in% root_terms), c("term_id", "name"), drop = FALSE]
if (!nrow(bp_terms)) stop("No non-root BP terms were found in terms.tsv")
term_exclusions <- unique(term_exclusions[, c("term_id", "reason")])
unknown_exclusions <- setdiff(term_exclusions$term_id, bp_terms$term_id)
if (length(unknown_exclusions)) {
  stop("term_exclusions.tsv contains IDs that are not non-root BP terms in terms.tsv: ",
       paste(unknown_exclusions, collapse = ", "))
}
bp_terms$excluded_by_frozen_protocol <- bp_terms$term_id %in% term_exclusions$term_id
memberships <- memberships[memberships$term_id %in% bp_terms$term_id, , drop = FALSE]
# Exclusions are part of the definition of a valid host-BP label.  Remove their
# memberships before defining U so units supported only by excluded labels do
# not silently remain in the testable denominator.
memberships <- memberships[!(memberships$term_id %in% term_exclusions$term_id), , drop = FALSE]

# U is the actual testable BP universe.  It is not the whole 41,139 nominal universe.
U_ids <- sort(unique(memberships$unit_id))
U <- universe[match(U_ids, universe$unit_id), , drop = FALSE]
if (!nrow(U)) stop("No eligible BP-annotated universe units after input filtering")

direction_values <- sort(unique(as.character(universe$direction[universe$opposite])))
if (!(cfg$direction_a %in% direction_values) || !(cfg$direction_b %in% direction_values)) {
  stop("Configured direction labels are not both present among opposite units. Available labels: ",
       paste(direction_values, collapse = " | "),
       ". Supply --direction-a and --direction-b explicitly.")
}

is_A_all <- universe$opposite & as.character(universe$direction) == cfg$direction_a
is_B_all <- universe$opposite & as.character(universe$direction) == cfg$direction_b
is_391_all <- universe$opposite
if (any(is_A_all & is_B_all)) stop("Direction A and B overlap")
if (sum(is_A_all) + sum(is_B_all) != sum(is_391_all)) {
  stop("Opposite units contain direction labels other than configured A/B; protocol requires exactly two classes.")
}
if (sum(is_391_all) != 391L || sum(is_A_all) != 348L || sum(is_B_all) != 43L ||
    sum(universe$opposite & universe$strict127) != 127L ||
    sum(is_A_all & universe$strict127) != 106L || sum(is_B_all & universe$strict127) != 21L) {
  stop("Input selection counts do not match the frozen protocol: expected 391/348/43 and strict 127/106/21. ",
       "Do not analyse until the selection ledger is reconciled.")
}

# Map selections onto U in a stable unit order.
sel391 <- U$opposite
selA <- U$opposite & as.character(U$direction) == cfg$direction_a
selB <- U$opposite & as.character(U$direction) == cfg$direction_b
strict391 <- U$opposite & U$strict127
strictA <- selA & U$strict127
strictB <- selB & U$strict127
if (sum(strictA) + sum(strictB) != sum(strict391)) stop("strict127 opposite units are not exclusively A/B")

term_sizes <- table(factor(memberships$term_id, levels = bp_terms$term_id))
term_sizes <- as.integer(term_sizes)
names(term_sizes) <- bp_terms$term_id
term_ok <- term_sizes >= 5L & term_sizes <= floor(0.95 * nrow(U)) &
  !(bp_terms$term_id %in% term_exclusions$term_id)
eligible_terms <- bp_terms$term_id[term_ok]
if (!length(eligible_terms)) stop("No BP terms satisfy the fixed size rule: 5 <= n <= 95% of U")
bp_terms$background_n <- term_sizes[bp_terms$term_id]
bp_terms$eligible <- bp_terms$term_id %in% eligible_terms

# Direct memberships as a unit x term logical incidence matrix.  This is a local
# implementation detail; exported counts always remain in sequence-unit units.
M <- matrix(FALSE, nrow = nrow(U), ncol = length(eligible_terms),
            dimnames = list(U$unit_id, eligible_terms))
idx_u <- match(memberships$unit_id, U$unit_id)
idx_t <- match(memberships$term_id, eligible_terms)
keep <- !is.na(idx_t) & !is.na(idx_u)
M[cbind(idx_u[keep], idx_t[keep])] <- TRUE

fisher_rows <- function(selected, alternative, analysis, q_is_decision = TRUE) {
  n_selected <- sum(selected)
  n_bg <- length(selected)
  if (n_selected == 0L || n_selected == n_bg) stop("Invalid selected-set size for ", analysis)
  ans <- vector("list", ncol(M))
  for (j in seq_len(ncol(M))) {
    a <- sum(selected & M[, j]); b <- n_selected - a
    c <- sum(!selected & M[, j]); d <- (n_bg - n_selected) - c
    tab <- matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE)
    # The requested p-value follows the analysis alternative.  Exact OR CIs are
    # always calculated two-sided, so A does not inherit a one-sided [low, Inf]
    # interval from its one-sided enrichment test.
    ft_p <- fisher.test(tab, alternative = alternative)
    ft_ci <- fisher.test(tab, alternative = "two.sided")
    ans[[j]] <- data.frame(
      analysis = analysis, term_id = colnames(M)[j], selected_n = n_selected,
      selected_term_n = a, background_n = n_bg, background_term_n = a + c,
      odds_ratio = unname(ft_ci$estimate), ci_low = ft_ci$conf.int[1], ci_high = ft_ci$conf.int[2],
      p_value = ft_p$p.value, stringsAsFactors = FALSE
    )
  }
  out <- do.call(rbind, ans)
  out$q_within_analysis <- p.adjust(out$p_value, method = "BH")
  out$decision_family <- q_is_decision
  out
}

# A = primary overall 391 ORA; B = composition contrast inside the 391 set.
res_A <- fisher_rows(sel391, "greater", "A_overall391_vs_annotated_background", TRUE)
if (sum(selB) >= 10L) {
  # Restrict the composition comparison to the selected 391 units; B is the
  # smaller directional class.  `selected` here encodes membership in class A.
  in391 <- sel391
  M_B <- M[in391, , drop = FALSE]
  selected_A_within391 <- selA[in391]
  old_M <- M; M <- M_B
  res_B <- fisher_rows(selected_A_within391, "two.sided", "B_directionA_vs_directionB_within391", TRUE)
  M <- old_M
} else {
  res_B <- NULL
}

# Combined BH is the only inferential q used for main-text decisions.  Separate
# within-analysis BH values are retained so the full multiplicity structure is visible.
primary <- if (is.null(res_B)) res_A else rbind(res_A, res_B)
primary$q_combined <- p.adjust(primary$p_value, method = "BH")
primary$enriched_class <- ifelse(primary$analysis == "A_overall391_vs_annotated_background", "overall391",
  ifelse(primary$odds_ratio >= 1, "directionA", "directionB"))
primary$enriched_class_term_n <- ifelse(primary$analysis == "A_overall391_vs_annotated_background",
  primary$selected_term_n,
  ifelse(primary$odds_ratio >= 1, primary$selected_term_n, primary$background_term_n - primary$selected_term_n))
primary$reportable_under_fixed_rule <- primary$q_combined < 0.05 &
  ((primary$analysis == "A_overall391_vs_annotated_background" & primary$odds_ratio > 1) |
   (primary$analysis == "B_directionA_vs_directionB_within391" & !is.na(primary$odds_ratio) & primary$odds_ratio != 1)) &
  primary$enriched_class_term_n >= 2L

# Nested strict-127 sensitivity outputs: all terms, no discovery calls.  Their q
# values are descriptive only and never used in reportable_under_fixed_rule.
sens_A127 <- fisher_rows(strict391, "greater", "S1_strict127_vs_annotated_background", FALSE)
if (sum(strictB) >= 10L) {
  in_strict <- strict391
  M_S <- M[in_strict, , drop = FALSE]
  old_M <- M; M <- M_S
  sens_B106v21 <- fisher_rows(strictA[in_strict], "two.sided", "S2_strictA_vs_strictB_within127", FALSE)
  M <- old_M
} else {
  sens_B106v21 <- NULL
}

# Quantile strata use rank-based quintiles, which remain defined under ties.
quintile <- function(x) {
  n <- length(x)
  as.integer(pmin(5L, floor(rank(x, ties.method = "average", na.last = "keep") / (n + 1) * 5) + 1L))
}
valid_bias <- is.finite(U$length) & U$length > 0 & is.finite(U$baseMean) & U$baseMean >= 0
if (!all(valid_bias)) {
  bad <- sum(!valid_bias)
  stop("Bias sensitivity not run: ", bad,
       " BP-annotated units lack a positive finite length or non-negative finite baseMean. Fix input rather than imputing.")
}
U$length_quintile <- quintile(U$length)
U$baseMean_quintile <- quintile(U$baseMean)
U$stratum <- interaction(U$length_quintile, U$baseMean_quintile, drop = TRUE, lex.order = TRUE)

make_term_by_stratum <- function(row_index) {
  # Rows=occupied strata, columns=eligible terms; each entry is a count of
  # annotated units carrying the term in that stratum.
  st <- droplevels(U$stratum[row_index])
  lv <- levels(st)
  z <- matrix(0L, nrow = length(lv), ncol = ncol(M), dimnames = list(lv, colnames(M)))
  for (s in seq_along(lv)) {
    ii <- which(st == lv[s])
    if (length(ii)) z[s, ] <- colSums(M[row_index[ii], , drop = FALSE])
  }
  z
}

set.seed(cfg$seed)
Bperm <- cfg$n_perm

# A permutation null: within each stratum retain exactly the observed number of
# selected 391 units.  For a term, the null count is a sum of hypergeometrics.
strata_all <- droplevels(U$stratum)
lev_all <- levels(strata_all)
N_all <- as.integer(table(strata_all)[lev_all])
n_Aall <- vapply(lev_all, function(s) sum(sel391 & strata_all == s), integer(1))
K_all <- make_term_by_stratum(seq_len(nrow(U)))
obs_A_count <- colSums(M[sel391, , drop = FALSE])
null_A <- matrix(0L, nrow = Bperm, ncol = ncol(M))
for (s in seq_along(lev_all)) {
  if (n_Aall[s] == 0L) next
  for (j in seq_len(ncol(M))) {
    k <- K_all[s, j]
    if (k == 0L || k == N_all[s]) {
      null_A[, j] <- null_A[, j] + if (k == N_all[s]) n_Aall[s] else 0L
    } else {
      null_A[, j] <- null_A[, j] + rhyper(Bperm, m = k, n = N_all[s] - k, k = n_Aall[s])
    }
  }
}
perm_A <- data.frame(term_id = colnames(M), observed_term_n = as.integer(obs_A_count),
                     perm_mean = colMeans(null_A), perm_sd = apply(null_A, 2L, sd),
                     p_empirical_enrichment = (1 + colSums(null_A >= rep(obs_A_count, each = Bperm))) / (Bperm + 1),
                     stringsAsFactors = FALSE)

# B permutation null: only within 391.  Each stratum preserves A and B class
# counts; the statistic is membership proportion(A) - membership proportion(B).
perm_B <- NULL
if (!is.null(res_B)) {
  in391 <- sel391
  strata_391 <- droplevels(U$stratum[in391])
  lev_391 <- levels(strata_391)
  N_391 <- as.integer(table(strata_391)[lev_391])
  a_391 <- vapply(lev_391, function(s) sum(selA[in391] & strata_391 == s), integer(1))
  K_391 <- make_term_by_stratum(which(in391))
  MA391 <- M[in391, , drop = FALSE]
  count_A_obs <- colSums(MA391[selA[in391], , drop = FALSE])
  count_total <- colSums(MA391)
  nA_total <- sum(selA); nB_total <- sum(selB)
  obs_diff <- count_A_obs / nA_total - (count_total - count_A_obs) / nB_total
  null_diff <- matrix(0, nrow = Bperm, ncol = ncol(M))
  # E0 is the analytic conditional expectation of the proportion difference
  # under the fixed stratum-specific A/B allocation.  It is generally not zero
  # when term membership and covariate strata are associated.
  expected_diff <- rep(0, ncol(M))
  for (s in seq_along(lev_391)) {
    for (j in seq_len(ncol(M))) {
      k <- K_391[s, j]
      e_x <- a_391[s] * k / N_391[s]
      expected_diff[j] <- expected_diff[j] + e_x / nA_total - (k - e_x) / nB_total
      if (a_391[s] == 0L) {
        # This stratum supplies no A unit but contributes the fixed B component
        # -k/nB.  It must not be skipped.
        x <- rep(0L, Bperm)
      } else if (a_391[s] == N_391[s]) {
        x <- rep(k, Bperm)
      } else {
        x <- rhyper(Bperm, m = k, n = N_391[s] - k, k = a_391[s])
      }
      # The total term count in the stratum is fixed, so B receives k-x.
      null_diff[, j] <- null_diff[, j] + x / nA_total - (k - x) / nB_total
    }
  }
  perm_B <- data.frame(term_id = colnames(M), observed_prop_difference_A_minus_B = obs_diff,
                       conditional_expectation = expected_diff,
                       observed_minus_expectation = obs_diff - expected_diff,
                       perm_mean = colMeans(null_diff), perm_sd = apply(null_diff, 2L, sd),
                       p_empirical_two_sided = (1 + colSums(abs(sweep(null_diff, 2L, expected_diff, "-")) >=
                                                            rep(abs(obs_diff - expected_diff), each = Bperm) - 1e-12)) / (Bperm + 1),
                       stringsAsFactors = FALSE)
}

# Attach terms and sensitivity labels, then write auditable outputs.
term_lookup <- bp_terms[, c("term_id", "name", "background_n", "eligible")]
names(term_lookup)[names(term_lookup) == "background_n"] <- "U_term_n"
primary <- merge(primary, term_lookup, by = "term_id", all.x = TRUE, sort = FALSE)
sens_A127 <- merge(sens_A127, term_lookup, by = "term_id", all.x = TRUE, sort = FALSE)
if (!is.null(sens_B106v21)) sens_B106v21 <- merge(sens_B106v21, term_lookup, by = "term_id", all.x = TRUE, sort = FALSE)
# Empirical p-values align with their own analysis only.  They receive an
# additional, separate BH correction over the same A+B 2T family; this does not
# replace Fisher q_combined, which remains the primary inferential decision.
primary$p_empirical <- NA_real_
primary$p_empirical[primary$analysis == "A_overall391_vs_annotated_background"] <-
  perm_A$p_empirical_enrichment[match(primary$term_id[primary$analysis == "A_overall391_vs_annotated_background"], perm_A$term_id)]
if (!is.null(perm_B)) {
  primary$p_empirical[primary$analysis == "B_directionA_vs_directionB_within391"] <-
    perm_B$p_empirical_two_sided[match(primary$term_id[primary$analysis == "B_directionA_vs_directionB_within391"], perm_B$term_id)]
}
primary$q_perm_combined <- NA_real_
has_perm <- is.finite(primary$p_empirical)
primary$q_perm_combined[has_perm] <- p.adjust(primary$p_empirical[has_perm], method = "BH")

coverage <- data.frame(
  metric = c("nominal_nonrisk_universe", "bp_annotated_testable_universe_U",
             "overall391_nominal", "overall391_in_U", "directionA_nominal", "directionA_in_U",
             "directionB_nominal", "directionB_in_U", "strict127_nominal", "strict127_in_U",
             "strictA_nominal", "strictA_in_U", "strictB_nominal", "strictB_in_U",
             "bp_terms_before_size_filter", "bp_terms_eligible_size_filter"),
  value = c(nrow(universe), nrow(U), sum(is_391_all), sum(sel391), sum(is_A_all), sum(selA),
            sum(is_B_all), sum(selB), sum(universe$opposite & universe$strict127), sum(strict391),
            sum(is_A_all & universe$strict127), sum(strictA), sum(is_B_all & universe$strict127), sum(strictB),
            nrow(bp_terms), length(eligible_terms)),
  stringsAsFactors = FALSE
)

strata_audit <- as.data.frame(table(U$stratum), stringsAsFactors = FALSE)
names(strata_audit) <- c("stratum", "U_n")
strata_audit$overall391_n <- vapply(strata_audit$stratum, function(s) sum(sel391 & U$stratum == s), integer(1))
strata_audit$directionA_n <- vapply(strata_audit$stratum, function(s) sum(selA & U$stratum == s), integer(1))
strata_audit$directionB_n <- vapply(strata_audit$stratum, function(s) sum(selB & U$stratum == s), integer(1))

write.table(coverage, file.path(cfg$output_dir, "00_annotation_coverage.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(term_lookup, file.path(cfg$output_dir, "01_bp_term_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(term_exclusions, file.path(cfg$output_dir, "01b_term_exclusions_applied.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(primary, file.path(cfg$output_dir, "02_primary_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sens_A127, file.path(cfg$output_dir, "03_sensitivity_strict127_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(sens_B106v21)) write.table(sens_B106v21, file.path(cfg$output_dir, "04_sensitivity_strictA_vs_strictB_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
if (!is.null(perm_B)) write.table(perm_B, file.path(cfg$output_dir, "05_permutation_B_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(perm_A, file.path(cfg$output_dir, "06_permutation_A_all_terms.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(strata_audit, file.path(cfg$output_dir, "07_bias_strata_audit.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

protocol <- list(
  script = "run_go_statistics.R",
  seed = cfg$seed, n_perm = Bperm,
  direction_A = cfg$direction_a, direction_B = cfg$direction_b,
  direct_BP_only = TRUE, go_ancestor_propagation = FALSE,
  frozen_term_exclusion_file = "term_exclusions.tsv",
  source_protein_GO_used = FALSE, candidate_domain_GO_used = FALSE,
  testable_background = "unique non-risk universe units with >=1 direct, non-root legacy host BP label",
  term_filter = "background membership >=5 and <=95% of U",
  primary_A = "one-sided Fisher: overall391 vs U",
  primary_B = if (is.null(res_B)) "not run: annotated direction B <10" else "two-sided Fisher: direction A vs B within annotated overall391",
  primary_q = "BH over all eligible A and B tests combined; q_combined is the sole main-text decision q",
  strict_sensitivity = "nested only; all-term output, no new discovery call",
  permutation_A = "term-wise conditional hypergeometric Monte Carlo within length quintile x baseMean quintile; preserve selected overall391 count per stratum",
  permutation_B = if (is.null(perm_B)) "not run" else "term-wise conditional hypergeometric Monte Carlo within annotated overall391 strata; exchange A/B labels preserving A count per stratum; two-sided deviations centered on analytic conditional expectation",
  permutation_q = "BH across all available A+B empirical term p-values (q_perm_combined); sensitivity only, separate from Fisher q_combined",
  caveat = "GO membership enrichment is not process activity, expression direction, or causal function"
)
json <- function(x) {
  # jsonlite is available in the approved R environment; use it only for a
  # machine-readable protocol manifest, not for statistical computation.
  jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required by the approved R environment for manifest output")
writeLines(json(protocol), file.path(cfg$output_dir, "00_protocol_manifest.json"))

cat("Completed frozen-protocol GO statistics. Outputs written to: ", normalizePath(cfg$output_dir), "\n", sep = "")

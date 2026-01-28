#' Run per-recruit SPR and YPR calculations
#'
#' Computes spawning potential ratio (SPR) and yield per recruit (YPR) across a
#' range of fishing mortality values for a single stock or multispecies complex
#' defined by \code{spr_input()}. The function returns species-specific and
#' combined (total) SPR(F) and YPR(F) curves, SPR-based fishing mortality
#' reference points (e.g., F40, F35), and optional decomposition diagnostics.
#'
#' @details
#' \strong{SPR and YPR.} For each species, unfished spawning biomass per recruit
#' (SBPR0) and fished SBPR(F) are computed from survivorship-at-age,
#' weight-at-age, and maturity-at-age. SPR is defined as \code{SPR(F) = SBPR(F)
#' / SBPR(0)}. YPR is computed using the Baranov catch equation in biomass units
#' and the age-specific fishery selectivity.
#'
#' \strong{Multispecies complexes.} For multispecies inputs, each species is
#' evaluated separately using its own life history schedules and recruitment
#' scaling (i.e., recruitment proportions in \code{spr_input()}). Combined SBPR
#' and yield are computed by summing across species.
#'
#' \strong{SPR reference points.} Fishing mortality reference points are reported for:
#' \itemize{
#'   \item the combined complex (\code{F_spr_total}), based on matching the combined SPR curve to
#'   each value in \code{spr_targets} using a nearest-grid approach, and
#'   \item each species individually (\code{F_spr_by_species}).
#' }
#'
#' \strong{Multispecies constraint.} When \code{multispecies_constraint =
#' "all_species"}, a constrained reference point
#' (\code{F_spr_total_constrained}) is calculated for each SPR target as the max
#' F in \code{F_vec} such that \emph{both} the combined SPR and each
#' species-specific SPR meet or exceed the target. If no value satisfies the
#' constraint, \code{NA} is returned for that target.
#'
#' \strong{Diagnostics and decomposition outputs.} When \code{diagnostics = TRUE}, the output includes:
#' \itemize{
#'   \item \code{diagnostics$plus_share}: the fraction of total SBPR contributed by the terminal age
#'   (often a plus group) at F = 0 and at the diagnostic reference F \code{Ftarget}, and
#'   \item \code{diagnostics$decomp}: an age-specific df containing survivorship, maturity,
#'   selectivity, weight-at-age, M-at-age, SBPR contributions, and SBPR removed
#'   (F0 minus Ftarget) for each species at F = 0 and F = Ftarget.
#' }
#' These diagnostics can be visualized with \code{plot_spr_decomp()}.
#'
#' @param x An object returned by \code{spr_input()} defining ages and species-specific life history schedules.
#' @param F_vec Numeric vector of fishing mortality values (F) to evaluate. Must be finite and >= 0.
#' @param spr_targets Numeric vector of SPR targets in (0,1), for example \code{c(0.40, 0.35)}.
#' @param multispecies_constraint Character string specifying whether to apply a multispecies constraint
#'   when selecting combined reference points. Options are:
#'   \describe{
#'     \item{\code{"none"}}{Unconstrained reference points are based on the combined SPR curve only.}
#'     \item{\code{"all_species"}}{Constrained combined reference points, all species must meet the target SPR.}
#'   }
#' @param diagnostics Logical. If TRUE, returns plus-group diagnostics and age-specific decompositions.
#'
#' @return A named list with the following elements:
#' \describe{
#'   \item{\code{F_vec}}{Fishing mortality grid used for calculations.}
#'   \item{\code{spr_targets}}{SPR targets requested.}
#'   \item{\code{multispecies_constraint}}{Constraint setting used.}
#'   \item{\code{SPR_total}}{Combined SPR(F) across species (length = \code{length(F_vec)}).}
#'   \item{\code{YPR_total}}{Combined YPR(F) across species (length = \code{length(F_vec)}).}
#'   \item{\code{SPR_by_species}}{Named list of species-specific SPR(F) vectors.}
#'   \item{\code{YPR_by_species}}{Named list of species-specific YPR(F) vectors.}
#'   \item{\code{F_spr_total}}{Named numeric vector of combined reference points (nearest-grid match),
#'     for example \code{F40_total} and \code{F35_total}.}
#'   \item{\code{F_spr_total_constrained}}{Named numeric vector of constrained combined reference points
#'     when \code{multispecies_constraint="all_species"}; otherwise \code{NULL}.}
#'   \item{\code{F_spr_by_species}}{Named list of species-specific reference points for each SPR target.}
#'   \item{\code{diagnostics}}{A list with elements:
#'     \itemize{
#'       \item \code{Ftarget}: the fishing mortality value used for diagnostics (defaults to the first target,
#'       and uses constrained \code{Ftarget} when available),
#'       \item \code{plus_share}: plus-group share table (or \code{NULL} if \code{diagnostics=FALSE}),
#'       \item \code{decomp}: age-specific decomposition data frame (or \code{NULL} if \code{diagnostics=FALSE}).
#'     }}
#'   }
#'
#' @examples
#' \dontrun{
#' # Assume inp is a spr_input() object (single species or multispecies).
#' # Run per-recruit calculations and request decomposition diagnostics.
#' spr_out <- run_spr(inp, diagnostics = TRUE,
#'                   multispecies_constraint = "all_species")
#'
#' # reference points
#' spr_out$F_spr_total
#' spr_out$F_spr_total_constrained
#' spr_out$F_spr_by_species
#'
#' # Plot SPR(F) curves, complex-level SPR and species-specific SPRs
#' # plot_spr_curves(spr_out, which = "total")
#' # plot_spr_curves(spr_out, which = "species")
#'
#' # Decomposition diagnostics
#' plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
#' plots$contrib
#' plots$removed
#' plots$survivorship
#'
#' # Plus-group diagnostic table (fraction of SBPR in terminal age)
#' spr_out$diagnostics$plus_share
#' }
#'
#' @export
run_spr <- function(
    x,
    F_vec = seq(0, 0.4, by = 0.001),
    spr_targets = c(0.40, 0.35),
    multispecies_constraint = c("none", "all_species"),
    diagnostics = TRUE
) {
  multispecies_constraint <- match.arg(multispecies_constraint)

  if (!is.list(x) || is.null(x$ages) || is.null(x$species)) stop("x must be an object returned by spr_input().")

  ages <- x$ages
  sp_names <- names(x$species)
  nsp <- length(sp_names)

  if (!is.numeric(F_vec) || any(!is.finite(F_vec)) || any(F_vec < 0)) stop("F_vec must be numeric, finite, and >= 0.")
  if (!is.numeric(spr_targets) || any(!is.finite(spr_targets)) || any(spr_targets <= 0) || any(spr_targets >= 1)) {
    stop("spr_targets must be numeric values in (0,1), e.g., c(0.40, 0.35).")
  }

  # solve for Fx such that SPR(Fx) = x over a defined F grid
  find_F_at_target <- function(SPR, target) {
    F_vec[which.min(abs(SPR - target))]
  }

  # get per recruit results for one species
  spr_engine_one <- function(sp, F) {
    M   <- sp$M_at_age
    m   <- sp$maturity_at_age
    sel <- sp$selex_at_age
    w   <- sp$wt_at_age
    R0  <- sp$R0
    n <- length(ages)

    # numbers at age
    N0 <- numeric(n)
    NF <- numeric(n)
    N0[1] <- R0
    NF[1] <- R0

    for (i in 2:n) {
      N0[i] <- N0[i - 1] * exp(-M[i - 1])
      NF[i] <- NF[i - 1] * exp(-(M[i - 1] + F * sel[i - 1]))
    }

    # plus-group
    if (isTRUE(x$use_plus_group)) {
      N0[n] <- N0[n - 1] * exp(-M[n - 1]) / (1 - exp(-M[n]))
      Z_last <- M[n] + F * sel[n]
      NF[n] <- NF[n - 1] * exp(-(M[n - 1] + F * sel[n - 1])) /
        (1 - exp(-Z_last))
    }

    # SPR components
    SB0 <- sum(N0 * w * m)
    SBF <- sum(NF * w * m)

    # Yield per recruit in wt
    Z <- M + F * sel
    Cn <- (F * sel / Z) * (1 - exp(-Z)) * NF
    Y  <- sum(Cn * w)

    # Per-recruit survivorship (scaled by R0)
    l0 <- N0 / R0
    lF <- NF / R0

    # Age-spec SBPR contributions
    contrib0 <- l0 * w * m
    contribF <- lF * w * m

    out <- list(
      SBPR0 = SB0,
      SBPR  = SBF,
      YPR   = Y,
      l0 = l0,
      lF = lF,
      contrib0 = contrib0,
      contribF = contribF,
      maturity = m,
      selectivity = sel,
      weight = w,
      M = M
    )

    if (isTRUE(x$use_plus_group)) {
      out$plus_share0 <- tail(contrib0, 1) / sum(contrib0)
      out$plus_shareF <- tail(contribF, 1) / sum(contribF)
    }

    out
  }

  # Run for each species
  by_sp <- stats::setNames(vector("list", nsp), sp_names)

  for (spn in sp_names) {
    sp <- x$species[[spn]]
    out_list <- lapply(F_vec, function(F) spr_engine_one(sp, F))

    SB0 <- out_list[[1]]$SBPR0
    SB  <- vapply(out_list, function(z) z$SBPR, numeric(1))
    Y   <- vapply(out_list, function(z) z$YPR,  numeric(1))
    SPR <- SB / SB0

    by_sp[[spn]] <- list(
      SPR = SPR,
      YPR = Y,
      SBPR0 = SB0,
      out_list = if (isTRUE(diagnostics)) out_list else NULL
    )
  }

  # Total SPR and YPR
  SB0_tot <- sum(vapply(sp_names, function(spn) by_sp[[spn]]$SBPR0, numeric(1)))
  SB_tot_byF <- rowSums(sapply(sp_names, function(spn) by_sp[[spn]]$SPR * by_sp[[spn]]$SBPR0))
  SPR_total  <- SB_tot_byF / SB0_tot
  YPR_total  <- rowSums(sapply(sp_names, function(spn) by_sp[[spn]]$YPR))

  # Unconstrained F targets (complex and by species)
  F_spr_total <- stats::setNames(
    vapply(spr_targets, function(t) find_F_at_target(SPR_total, t), numeric(1)),
    paste0("F", round(100 * spr_targets), "_total")
  )

  F_spr_by_species <- lapply(sp_names, function(spn) {
    stats::setNames(
      vapply(spr_targets, function(t) find_F_at_target(by_sp[[spn]]$SPR, t), numeric(1)),
      paste0("F", round(100 * spr_targets), "_", spn)
    )
  })
  names(F_spr_by_species) <- sp_names

  # Optional constrained F targets
  F_spr_total_constrained <- NULL
  if (multispecies_constraint == "all_species") {
    F_spr_total_constrained <- stats::setNames(
      rep(NA_real_, length(spr_targets)),
      paste0("F", round(100 * spr_targets), "_total_constrained")
    )

    SPR_mat <- sapply(sp_names, function(spn) by_sp[[spn]]$SPR)

    for (i in seq_along(spr_targets)) {
      t <- spr_targets[i]
      ok <- (SPR_total >= t) & apply(SPR_mat >= t, 1, all)
      F_spr_total_constrained[i] <- if (!any(ok)) NA_real_ else max(F_vec[ok])
    }
  }

  # define default diagnostic F (uses constrained if available, else unconstrained)
  Fdiag <- unname(F_spr_total[[1]])
  if (multispecies_constraint == "all_species" &&
      !is.null(F_spr_total_constrained) &&
      is.finite(F_spr_total_constrained[[1]])) {
    Fdiag <- unname(F_spr_total_constrained[[1]])
  }

  plus_share <- NULL
  decomp <- NULL

  if (isTRUE(diagnostics)) {

    if (isTRUE(x$use_plus_group)) {
      plus_share <- do.call(rbind, lapply(sp_names, function(spn) {
        sp_out0 <- by_sp[[spn]]$out_list[[1]]
        sp_outF <- by_sp[[spn]]$out_list[[which.min(abs(F_vec - Fdiag))]]
        data.frame(
          species = spn,
          Flevel = c("F0", "Ftarget"),
          F = c(0, Fdiag),
          plus_share = c(sp_out0$plus_share0, sp_outF$plus_shareF),
          stringsAsFactors = FALSE
        )
      }))
    }

    decomp <- do.call(rbind, lapply(sp_names, function(spn) {
      sp_out0 <- by_sp[[spn]]$out_list[[1]]
      sp_outF <- by_sp[[spn]]$out_list[[which.min(abs(F_vec - Fdiag))]]

      rbind(
        data.frame(
          species = spn,
          age = ages,
          Flevel = "F0",
          F = 0,
          survivorship = sp_out0$l0,
          maturity = sp_out0$maturity,
          selectivity = sp_out0$selectivity,
          weight = sp_out0$weight,
          M = sp_out0$M,
          contrib = sp_out0$contrib0,
          stringsAsFactors = FALSE
        ),
        data.frame(
          species = spn,
          age = ages,
          Flevel = "Ftarget",
          F = Fdiag,
          survivorship = sp_outF$lF,
          maturity = sp_outF$maturity,
          selectivity = sp_outF$selectivity,
          weight = sp_outF$weight,
          M = sp_outF$M,
          contrib = sp_outF$contribF,
          stringsAsFactors = FALSE
        )
      )
    }))

    decomp_rem <- stats::reshape(
      decomp[, c("species","age","Flevel","contrib")],
      timevar = "Flevel",
      idvar = c("species","age"),
      direction = "wide"
    )
    decomp_rem$removed <- decomp_rem$contrib.F0 - decomp_rem$contrib.Ftarget

    decomp <- merge(
      decomp,
      decomp_rem[, c("species","age","removed")],
      by = c("species","age"),
      all.x = TRUE
    )
  }

  out <- list(
    F_vec = F_vec,
    spr_targets = spr_targets,
    multispecies_constraint = multispecies_constraint,
    SPR_total = SPR_total,
    YPR_total = YPR_total,
    SPR_by_species = lapply(by_sp, `[[`, "SPR"),
    YPR_by_species = lapply(by_sp, `[[`, "YPR"),
    F_spr_total = F_spr_total,
    F_spr_total_constrained = F_spr_total_constrained,
    F_spr_by_species = F_spr_by_species,
    diagnostics = list(
      Ftarget = Fdiag,
      plus_share = plus_share,
      decomp = decomp
    )
  )

  # order by species
  attr(out, "species") <- sp_names
  out
}


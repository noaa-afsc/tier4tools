#' Sensitivity analysis for SPR reference points
#'
#' Runs a grid-based sensitivity analysis for Tier 4 SPR calculations by varying
#' a user-defined set of influential inputs. Available options include: scalar
#' natural mortality (M), logistic maturity \eqn{a_{50}}, knife-edge selectivity
#' \eqn{a_c}, and (for 2-species complexes) recruitment proportions.
#'
#' @details This function modifies a base \code{spr_input()} object
#'   \code{x_base} by updating one or more schedules and re-running
#'   \code{run_spr()} for each scenario. The result is returned as a tidy data
#'   frame suitable for sensitivity plots or tables.
#'
#' \strong{Important limitation:} All sensitivity sweep dimensions are applied
#' to \emph{all species simultaneously}. For example, sweeping
#' \code{maturity_a50} sets the maturity ogive midpoint to the same value for
#' every species in the input, and sweeping \code{selex_ac} sets the same
#' knife-edge \eqn{a_c} for every species. Species-specific sensitivity sweeps
#' (varying a parameter for one species while holding others fixed) are not
#' currently implemented.
#'
#' \strong{What can be varied:}
#' \itemize{
#'   \item \code{M}: scalar \eqn{M} applied at all ages (replaces existing \code{M_at_age}).
#'   \item \code{maturity_a50}: maturity ogive midpoint \eqn{a_{50}} for a logistic maturity schedule.
#'   \item \code{selex_ac}: knife-edge selectivity age \eqn{a_c}.
#'   \item \code{prop}: recruitment proportion assigned to the first species (2-species inputs only).
#' }
#'
#' \strong{Outputs:} The returned df includes combined (total) and (optionally)
#' constrained SPR reference points (e.g., F40, F35), total YPR evaluated at
#' F40, and plus group diagnostics as separate columns for each species.
#'
#' @param x_base A \code{spr_input()} object used as the base scenario.
#' @param M Numeric vector of scalar M values to apply. If \code{NULL}, M is not
#'   varied.
#' @param maturity_a50 Numeric vector of logistic maturity a50 values. If
#'   \code{NULL}, not varied.
#' @param selex_ac Numeric vector of knife-edge selectivity a_c values. If
#'   \code{NULL}, not varied.
#' @param prop Numeric vector of recruitment proportions for 2-species inputs.
#'   Interpreted as the proportion assigned to the first species in
#'   \code{names(x_base$species)}. Ignored for 1 species.
#' @param F_vec Fishing mortality grid passed to \code{run_spr()}.
#' @param spr_targets SPR targets passed to \code{run_spr()} (default
#'   \code{c(0.40, 0.35)}).
#' @param multispecies_constraint Passed to \code{run_spr()} (\code{"none"} or
#'   \code{"all_species"}).
#' @param use_constrained_for_yield Logical. If TRUE and constrained F40 exists,
#'   evaluate \code{YPR_total} at constrained F40; otherwise evaluate at
#'   unconstrained F40.
#' @param digits Digits used for formatting \code{scenario_label}.
#'
#' @return A df with one row per scenario. Columns include:
#' \describe{
#'   \item{\code{scenario_id}}{Unique integer id.}
#'   \item{\code{M}, \code{maturity_a50}, \code{selex_ac}, \code{prop}}{Scenario inputs (NA means base).}
#'   \item{\code{F40_total}, \code{F35_total}}{Combined reference points (unconstrained).}
#'   \item{\code{F40_total_constrained}, \code{F35_total_constrained}}{Combined constrained reference points
#'     when \code{multispecies_constraint="all_species"} (otherwise NA).}
#'   \item{\code{YPR_at_F40}}{Total yield per recruit evaluated at F40 (constrained if requested).}
#'   \item{\code{plus_share_<sp>_F0}}{Plus-group share at F=0 for each species \code{<sp>}.}
#'   \item{\code{plus_share_<sp>_Ftarget}}{Plus-group share at the diagnostic Ftarget for each species.}
#'   \item{\code{scenario_label}}{Readable label concatenating scenario inputs.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- spr_sensitivity(
#'   x_base,
#'   M = seq(0.03, 0.06, by = 0.005),
#'   maturity_a50 = seq(10, 20, by = 2),
#'   selex_ac = c(8, 10, 12),
#'   prop = seq(0.3, 0.7, by = 0.1),
#'   F_vec = seq(0, 0.4, by = 0.002),
#'   multispecies_constraint = "all_species"
#' )
#'
#' # Heatmap: F40_total across (M, a50)
#' plot_sensitivity_heatmap(out, x = "M", y = "maturity_a50", fill = "F40_total")
#' }
#'
#' @export
spr_sensitivity <- function(
    x_base,
    M = NULL,
    maturity_a50 = NULL,
    selex_ac = NULL,
    prop = NULL,
    F_vec = seq(0, 0.4, by = 0.002),
    spr_targets = c(0.40, 0.35),
    multispecies_constraint = c("none", "all_species"),
    use_constrained_for_yield = TRUE,
    digits = 3
) {
  multispecies_constraint <- match.arg(multispecies_constraint)

  if (!is.list(x_base) || is.null(x_base$ages) || is.null(x_base$species)) {
    stop("x_base must be an object returned by spr_input().")
  }
  sp_names <- names(x_base$species)
  nsp <- length(sp_names)
  ages <- x_base$ages
  nages <- length(ages)

  any_sweep <- !is.null(M) || !is.null(maturity_a50) || !is.null(selex_ac) || !is.null(prop)
  if (nsp > 1 && any_sweep) {
    message(
      "spr_sensitivity(): Sweep values are applied to all species simultaneously. ",
      "Species-specific sweeps are not currently supported."
    )
  }

  # Baseline values
  base_M_scalar <- NA_real_
  {
    m_by_sp <- lapply(sp_names, function(spn) as.numeric(x_base$species[[spn]]$M_at_age))
    const_by_sp <- vapply(m_by_sp, function(v) {
      u <- unique(v[is.finite(v)])
      length(u) == 1
    }, logical(1))
    if (all(const_by_sp)) {
      vals <- vapply(m_by_sp, function(v) unique(v[is.finite(v)])[1], numeric(1))
      if (length(unique(round(vals, 12))) == 1) base_M_scalar <- as.numeric(vals[1])
    }
  }

  base_maturity_a50 <- NA_real_
  {
    a50s <- vapply(sp_names, function(spn) {
      ms <- x_base$species[[spn]]$maturity_spec
      if (is.list(ms) && identical(ms$type, "logistic") && is.finite(as.numeric(ms$a50))) {
        as.numeric(ms$a50)
      } else {
        NA_real_
      }
    }, numeric(1))
    if (all(is.finite(a50s)) && length(unique(round(a50s, 12))) == 1) base_maturity_a50 <- a50s[1]
  }

  base_selex_ac <- NA_real_
  {
    acs <- vapply(sp_names, function(spn) {
      ss <- x_base$species[[spn]]$spec$selectivity
      if (is.list(ss) && identical(ss$type, "knife_edge") && is.finite(as.numeric(ss$a_c))) {
        as.numeric(ss$a_c)
      } else {
        NA_real_
      }
    }, numeric(1))
    if (all(is.finite(acs)) && length(unique(round(acs, 12))) == 1) base_selex_ac <- acs[1]
  }

  base_prop <- NA_real_
  if (nsp == 2 && !is.null(x_base$rec_prop) && !is.null(names(x_base$rec_prop))) {
    if (sp_names[1] %in% names(x_base$rec_prop)) base_prop <- as.numeric(x_base$rec_prop[[sp_names[1]]])
  }

  # Sweep vectors
  M_vec    <- if (is.null(M)) NA_real_ else unique(c(NA_real_, as.numeric(M)))
  a50_vec  <- if (is.null(maturity_a50)) NA_real_ else unique(c(NA_real_, as.numeric(maturity_a50)))
  ac_vec   <- if (is.null(selex_ac)) NA_real_ else unique(c(NA_real_, as.numeric(selex_ac)))
  prop_vec <- if (is.null(prop)) NA_real_ else unique(c(NA_real_, as.numeric(prop)))

  if (!is.null(M) && any(!is.finite(M_vec[is.finite(M_vec)]) | M_vec[is.finite(M_vec)] <= 0)) stop("M must be finite and > 0.")
  if (!is.null(maturity_a50) && any(!is.finite(a50_vec[is.finite(a50_vec)]))) stop("maturity_a50 must be finite.")
  if (!is.null(selex_ac) && any(!is.finite(ac_vec[is.finite(ac_vec)]))) stop("selex_ac must be finite.")
  if (!is.null(prop)) {
    if (nsp != 2) stop("prop sweep is currently supported only for 2-species inputs.")
    if (any(!is.finite(prop_vec[is.finite(prop_vec)]) | prop_vec[is.finite(prop_vec)] < 0 | prop_vec[is.finite(prop_vec)] > 1)) {
      stop("prop must be in [0,1].")
    }
    for (spn in sp_names) {
      if (is.null(x_base$species[[spn]]$R0_base) || !is.finite(as.numeric(x_base$species[[spn]]$R0_base))) {
        stop("To sweep prop, store baseline recruitment scaling as x$species[[sp]]$R0_base in spr_input().")
      }
    }
  }
  if (!is.null(maturity_a50)) {
    for (spn in sp_names) {
      if (is.null(x_base$species[[spn]]$maturity_delta) || !is.finite(x_base$species[[spn]]$maturity_delta)) {
        stop("To sweep maturity_a50, each species must store a finite maturity_delta (x$species[[sp]]$maturity_delta).")
      }
    }
  }

  grid <- expand.grid(
    M = M_vec,
    maturity_a50 = a50_vec,
    selex_ac = ac_vec,
    prop = prop_vec,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  if (nrow(grid) == 0) stop("No scenarios generated.")

  apply_scenario <- function(x, row) {
    if (is.finite(row$M)) {
      for (spn in sp_names) x$species[[spn]]$M_at_age <- rep(row$M, nages)
    }
    if (is.finite(row$maturity_a50)) {
      for (spn in sp_names) {
        delta <- x$species[[spn]]$maturity_delta
        a50 <- row$maturity_a50
        m <- 1 / (1 + exp((ages - a50) / delta))
        x$species[[spn]]$maturity_at_age <- pmin(1, pmax(0, as.numeric(m)))
      }
    }
    if (is.finite(row$selex_ac)) {
      a_c <- row$selex_ac
      for (spn in sp_names) x$species[[spn]]$selex_at_age <- as.numeric(ifelse(ages >= a_c, 1, 0))
    }
    if (is.finite(row$prop)) {
      p <- row$prop
      x$species[[sp_names[1]]]$R0 <- x$species[[sp_names[1]]]$R0_base * p
      x$species[[sp_names[2]]]$R0 <- x$species[[sp_names[2]]]$R0_base * (1 - p)
      x$rec_prop <- c(stats::setNames(p, sp_names[1]), stats::setNames(1 - p, sp_names[2]))
      x$rec_prop <- x$rec_prop[sp_names]
    }
    x
  }

  # returns numeric value from named object
  get_named <- function(x, nm) {
    if (is.null(x) || is.null(names(x)) || !(nm %in% names(x))) return(NA_real_)
    as.numeric(x[[nm]])
  }

  # returns the sbpr contribution for the plus grp
  get_plus_share <- function(ps, spn, flevel) {
    if (is.null(ps) || !is.data.frame(ps) || nrow(ps) == 0) return(NA_real_)
    v <- ps$plus_share[ps$species == spn & ps$Flevel == flevel]
    if (length(v) == 0) NA_real_ else as.numeric(v[1])
  }

  # species-specific F columns from run_spr()
  flatten_F_by_species <- function(F_list) {
    if (is.null(F_list) || !is.list(F_list)) return(list())
    out <- list()
    for (spn in names(F_list)) {
      v <- F_list[[spn]]
      if (is.null(names(v))) next
      for (nm in names(v)) out[[nm]] <- as.numeric(v[[nm]])
    }
    out
  }

  # used in calc_yield_share() below - uses linear approximation to get YPR at F target F*
  .get_at_F <- function(F_vec, y_vec, F_value, method = c("linear", "nearest")) {
    method <- match.arg(method)
    if (!is.finite(F_value)) return(NA_real_)
    if (length(F_vec) != length(y_vec)) stop("F_vec and y_vec must have same length.")
    if (method == "nearest") return(y_vec[which.min(abs(F_vec - F_value))])
    stats::approx(x = F_vec, y = y_vec, xout = F_value, rule = 2)$y
  }
  # yield shares at a given F* using x$species[[sp]]$R0 as weights
  calc_yield_share_at_F <- function(spr_out, xi, F_star) {
    if (!is.finite(F_star)) {
      y_sp <- stats::setNames(rep(NA_real_, length(sp_names)), sp_names)
      return(list(YPR_sp = y_sp, YPR_total = NA_real_, share_sp = y_sp))
    }

    # evaluate species YPR at F*
    y_raw <- vapply(sp_names, function(spn) {
      .get_at_F(spr_out$F_vec, spr_out$YPR_by_species[[spn]], F_star, method = "linear")
    }, numeric(1))

    # recruitment weights from R0 (proportional)
    R0 <- vapply(sp_names, function(spn) as.numeric(xi$species[[spn]]$R0), numeric(1))
    if (any(!is.finite(R0)) || sum(R0) <= 0) R0 <- rep(1, length(R0))
    w <- R0 / sum(R0)

    y_w <- y_raw * w
    y_tot <- sum(y_w, na.rm = TRUE)
    share <- if (is.finite(y_tot) && y_tot > 0) y_w / y_tot else rep(NA_real_, length(y_w))

    list(
      YPR_sp = stats::setNames(y_w, sp_names),
      YPR_total = y_tot,
      share_sp = stats::setNames(share, sp_names)
    )
  }

  # Run scenarios
  res <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    xi <- apply_scenario(x_base, grid[i, ])

    spr_out <- run_spr(
      xi,
      F_vec = F_vec,
      spr_targets = spr_targets,
      multispecies_constraint = multispecies_constraint,
      diagnostics = TRUE
    )

    # total F targets (unconstrained)
    nm40 <- paste0("F", round(100 * spr_targets[1]), "_total")
    nm35 <- paste0("F", round(100 * spr_targets[2]), "_total")
    F40_total <- get_named(spr_out$F_spr_total, nm40)
    F35_total <- if (length(spr_targets) >= 2) get_named(spr_out$F_spr_total, nm35) else NA_real_

    # constrained total F targets if computed
    nm40c <- paste0("F", round(100 * spr_targets[1]), "_total_constrained")
    nm35c <- paste0("F", round(100 * spr_targets[2]), "_total_constrained")
    F40_con <- get_named(spr_out$F_spr_total_constrained, nm40c)
    F35_con <- if (length(spr_targets) >= 2) get_named(spr_out$F_spr_total_constrained, nm35c) else NA_real_

    # YPR_at_F40
    F_for_y <- F40_total
    if (isTRUE(use_constrained_for_yield) && is.finite(F40_con)) F_for_y <- F40_con
    idx_y <- which.min(abs(spr_out$F_vec - F_for_y))
    YPR_at_F40 <- spr_out$YPR_total[idx_y]

    # plus shares by species
    ps <- spr_out$diagnostics$plus_share
    plus_cols <- list()
    for (spn in sp_names) {
      plus_cols[[paste0("plus_share_", spn, "_F0")]] <- get_plus_share(ps, spn, "F0")
      plus_cols[[paste0("plus_share_", spn, "_Ftarget")]] <- get_plus_share(ps, spn, "Ftarget")
    }

    # species-specific ref point columns (e.g., F40_RE, F35_BS)
    F_by_sp_cols <- flatten_F_by_species(spr_out$F_spr_by_species)

    # yield shares at ref pts for each target (unconstrained + constrained)
    ycols <- list()
    for (t in spr_targets) {
      t_int <- round(100 * t)

      # unconstrained total reference F
      nm_tot <- paste0("F", t_int, "_total")
      F_star <- get_named(spr_out$F_spr_total, nm_tot)
      ys <- calc_yield_share_at_F(spr_out, xi, F_star)

      for (spn in sp_names) {
        ycols[[paste0("YPR", t_int, "_", spn)]] <- as.numeric(ys$YPR_sp[[spn]])
        ycols[[paste0("share", t_int, "_", spn)]] <- as.numeric(ys$share_sp[[spn]])
      }
      ycols[[paste0("YPR", t_int, "_total")]] <- as.numeric(ys$YPR_total)

      # constrained total reference F (if available)
      nm_totc <- paste0("F", t_int, "_total_constrained")
      F_star_c <- get_named(spr_out$F_spr_total_constrained, nm_totc)
      ys_c <- calc_yield_share_at_F(spr_out, xi, F_star_c)

      for (spn in sp_names) {
        ycols[[paste0("YPR", t_int, "_", spn, "_constrained")]] <- as.numeric(ys_c$YPR_sp[[spn]])
        ycols[[paste0("share", t_int, "_", spn, "_constrained")]] <- as.numeric(ys_c$share_sp[[spn]])
      }
      ycols[[paste0("YPR", t_int, "_total_constrained")]] <- as.numeric(ys_c$YPR_total)
    }

    base_row <- data.frame(
      scenario_id = i,
      M = if (is.finite(grid$M[i])) grid$M[i] else NA_real_,
      maturity_a50 = if (is.finite(grid$maturity_a50[i])) grid$maturity_a50[i] else NA_real_,
      selex_ac = if (is.finite(grid$selex_ac[i])) grid$selex_ac[i] else NA_real_,
      prop = if (is.finite(grid$prop[i])) grid$prop[i] else NA_real_,
      F40_total = F40_total,
      F35_total = F35_total,
      F40_total_constrained = F40_con,
      F35_total_constrained = F35_con,
      YPR_at_F40 = YPR_at_F40,
      stringsAsFactors = FALSE
    )

    res[[i]] <- cbind(
      base_row,
      as.data.frame(plus_cols, check.names = FALSE),
      as.data.frame(F_by_sp_cols, check.names = FALSE),
      as.data.frame(ycols, check.names = FALSE)
    )
  }

  out <- do.call(rbind, res)

  # format scenario label text
  fmt <- function(z) ifelse(is.na(z), "base", formatC(z, format = "f", digits = digits))
  prop_txt <- if (nsp == 2 && !is.null(prop)) paste0(", prop_", sp_names[1], "=", fmt(out$prop)) else ""
  out$scenario_label <- paste0(
    "M=", fmt(out$M),
    ", a50=", fmt(out$maturity_a50),
    ", ac=", fmt(out$selex_ac),
    prop_txt
  )

  out$M_value <- ifelse(is.na(out$M), base_M_scalar, out$M)
  out$maturity_a50_value <- ifelse(is.na(out$maturity_a50), base_maturity_a50, out$maturity_a50)
  out$selex_ac_value <- ifelse(is.na(out$selex_ac), base_selex_ac, out$selex_ac)
  out$prop_value <- ifelse(is.na(out$prop), base_prop, out$prop)

  attr(out, "base_values") <- list(
    M = base_M_scalar,
    maturity_a50 = base_maturity_a50,
    selex_ac = base_selex_ac,
    prop = base_prop
  )
  attr(out, "species") <- sp_names

  out
}

#' Plot sensitivity heatmap for SPR outputs
#'
#' Heatmap based on \code{spr_sensitivity()} output for a two-parameter sweep
#' with the base values overlayed. Uses viridis by default when available.
#'
#' @details
#' This function is intended for visualization of structured parameter sweeps
#' where \code{x} and \code{y} are numeric columns in the sensitivity table and
#' \code{fill} is a numeric output (for example \code{"F40_total"},
#' \code{"F40_total_constrained"}, \code{"YPR_at_F40"}, or
#' \code{"plus_share_RE_Ftarget"}).
#'
#' @param sens_df Data frame returned by \code{spr_sensitivity()}.
#' @param x Name of the column to use for the x-axis (character).
#' @param y Name of the column to use for the y-axis (character).
#' @param fill Name of the column to use for tile fill (character).
#' @param facet Optional column name to facet by (character).
#' @param title Optional plot title.
#' @param fill_scale Character, either \code{"base"} or \code{"viridis"} (default).
#' @param base_marker Logical. If TRUE, overlays a base-case marker when base values are available.
#' @param base_point_size Numeric, size of base case mark (X).
#' @param base_point_color Outline color for base case mark, defaults to a colorblind-friendly orange.
#' @param base_point_stroke Stroke width for base case mark.
#' @param base_label Logical. If TRUE, labels the base marker with "Base".
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # Heatmap of unconstrained F40 across (M, a50)
#' plot_sensitivity_heatmap(out, x = "M", y = "maturity_a50", fill = "F40_total")
#' }
#'
#' @export
plot_sensitivity_heatmap <- function(
    sens_df,
    x,
    y,
    fill,
    facet = NULL,
    title = NULL,
    fill_scale = c("viridis", "base"),
    base_marker = TRUE,
    base_point_size = 5,
    base_point_color = "#D55E00",
    base_point_stroke = 1.4,
    base_label = FALSE
) {
  fill_scale <- match.arg(fill_scale)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
  if (!requireNamespace("rlang", quietly = TRUE)) stop("Please install rlang.")

  if (!is.data.frame(sens_df)) stop("sens_df must be a data.frame from spr_sensitivity().")
  for (nm in c(x, y, fill)) {
    if (!nm %in% names(sens_df)) stop("Column not found in sens_df: ", nm)
  }
  if (!is.null(facet) && !facet %in% names(sens_df)) stop("Facet column not found in sens_df: ", facet)

  df <- sens_df

  x_plot <- if (paste0(x, "_value") %in% names(df)) paste0(x, "_value") else x
  y_plot <- if (paste0(y, "_value") %in% names(df)) paste0(y, "_value") else y

  # Keep only rows that can be drawn as tiles (x/y must exist)
  keep <- is.finite(as.numeric(df[[x_plot]])) & is.finite(as.numeric(df[[y_plot]]))
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) stop("No finite rows to plot after filtering x/y.")

  # Coerce x/y to ordered factors so tiles plot nicely
  make_ordered_factor <- function(v) {
    u <- sort(unique(as.numeric(v[is.finite(as.numeric(v))])))
    factor(as.numeric(v), levels = u, ordered = TRUE)
  }
  df[[x_plot]] <- make_ordered_factor(df[[x_plot]])
  df[[y_plot]] <- make_ordered_factor(df[[y_plot]])

  xq <- rlang::sym(x_plot)
  yq <- rlang::sym(y_plot)
  fq <- rlang::sym(fill)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = !!xq, y = !!yq, fill = !!fq)) +
    ggplot2::geom_tile() +
    ggplot2::labs(
      x = x,
      y = y,
      fill = fill,
      title = title
    ) +
    ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(mult = 0, add = 0))

  if (!is.null(facet)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet)))
  }

  # Fill scale
  if (fill_scale == "viridis") {
    if (requireNamespace("viridisLite", quietly = TRUE)) {
      p <- p + ggplot2::scale_fill_viridis_c(na.value = "transparent")
    } else {
      message("plot_sensitivity_heatmap(): 'viridisLite' not installed, using ggplot2 default fill scale.")
    }
  }

  # Overlay base case values
  if (isTRUE(base_marker)) {

    # Get base values from attribute of spr_sensitivity()
    bv <- attr(sens_df, "base_values", exact = TRUE)

    # Otherwise infer base vals from the row where all sweep inputs are NA
    sweep_cols <- intersect(c("M", "maturity_a50", "selex_ac", "prop"), names(sens_df))
    is_base_row <- rep(FALSE, nrow(sens_df))
    if (length(sweep_cols) > 0) {
      is_base_row <- apply(sens_df[, sweep_cols, drop = FALSE], 1, function(z) all(is.na(z)))
    }

    get_base_numeric <- function(param) {
      # attribute base_values
      if (is.list(bv) && !is.null(bv[[param]]) && is.finite(as.numeric(bv[[param]]))) {
        return(as.numeric(bv[[param]]))
      }
      # base row in the table, use *_value if it exists
      if (any(is_base_row)) {
        base_idx <- which(is_base_row)[1]
        cand <- paste0(param, "_value")
        if (cand %in% names(sens_df)) {
          v <- sens_df[[cand]][base_idx]
        } else if (param %in% names(sens_df)) {
          v <- sens_df[[param]][base_idx]
        } else {
          v <- NA_real_
        }
        if (is.finite(as.numeric(v))) return(as.numeric(v))
      }
      NA_real_
    }

    bx <- get_base_numeric(x)
    by <- get_base_numeric(y)

    if (is.finite(bx) && is.finite(by)) {
      # Convert numeric base to factor levels used in plotting
      bx_f <- factor(bx, levels = levels(df[[x_plot]]), ordered = TRUE)
      by_f <- factor(by, levels = levels(df[[y_plot]]), ordered = TRUE)

      base_df <- data.frame(.bx = bx_f, .by = by_f)

      # If faceting
      if (!is.null(facet) && facet %in% names(df)) {
        base_df[[facet]] <- df[[facet]][1]
      }

      # X for base
      p <- p +
        ggplot2::geom_point(
          data = base_df,
          ggplot2::aes(x = .bx, y = .by),
          inherit.aes = FALSE,
          shape = 4,
          size = base_point_size,
          stroke = base_point_stroke,
          color = base_point_color
        )

      if (isTRUE(base_label)) {
        p <- p +
          ggplot2::geom_text(
            data = base_df,
            ggplot2::aes(x = .bx, y = .by, label = "Base"),
            inherit.aes = FALSE,
            vjust = -0.8,
            color = base_point_color
          )
      }
    } else {
      message("plot_sensitivity_heatmap(): Base values not available, base marker not plotted.")
    }
  }

  p + theme_tier4()
}


#' Plot multispecies sensitivity to recruitment proportions
#'
#' Creates report-ready plots summarizing how multispecies (complex-level) SPR reference points
#' change across recruitment proportion assumptions for two-species complexes. The function is
#' designed to work directly with the output of \code{spr_sensitivity()} run under two settings:
#' (1) an unconstrained run where complex reference points are derived from the combined SPR curve
#' (\code{multispecies_constraint = "none"}), and (2) a constrained run where all species must meet
#' the target SPR (\code{multispecies_constraint = "all_species"}).
#'
#' @details
#' \strong{Default behavior.} The main plot is a single panel where:
#' \itemize{
#'   \item color = SPR target (for example, 0.40 vs 0.35), and
#'   \item linetype = the reference point method (unconstrained vs constrained).
#' }
#'
#' \strong{Delta option.} When \code{show_delta = TRUE}, the unconstrained values are plotted as a 0 baseline and the constrained values
#' shows the constraint-induced change at each recruitment proportion and target:
#' \deqn{\Delta F = F_{\mathrm{method}} - F_{\mathrm{unconstrained}}.}
#'
#' \strong{Optional diagnostics.}
#' \itemize{
#'   \item When \code{include_species = TRUE}, the function also attempts to return a limiting-species
#'   panel based on species-specific reference point columns produced by \code{spr_sensitivity()}
#'   (for example \code{F40_RE}, \code{F40_BS}). The limiting species is defined as the species with
#'   the smallest species-specific reference point at each recruitment proportion and target.
#'   \item When \code{include_yield_share = TRUE}, the function attempts to return an implied yield-share
#'   panel using yield-share columns added by \code{spr_sensitivity()}, for example
#'   \code{share40_RE}, \code{share40_BS}, and constrained equivalents
#'   \code{share40_RE_constrained}, \code{share40_BS_constrained}. Yield shares are interpreted as
#'   recruitment-weighted YPR contributions at the reference point F.
#' }
#' If the required columns are not present, the corresponding optional plots are returned as \code{NULL}.
#'
#' @param sens_uncon A data frame returned by \code{spr_sensitivity()} with
#'   \code{multispecies_constraint = "none"}.
#' @param sens_con A data frame returned by \code{spr_sensitivity()} with
#'   \code{multispecies_constraint = "all_species"}.
#' @param targets Numeric vector of SPR targets to plot (values in (0, 1)),
#'   for example \code{c(0.40, 0.35)}.
#' @param prop_col Character. Name of the recruitment proportion column to use on the x-axis.
#'   Defaults to \code{"prop_value"} (recommended), but can be set to \code{"prop"} if desired.
#' @param show_delta Logical. If FALSE (default), plots absolute reference point F.
#'   If TRUE, plots deviations from the unconstrained reference point (\eqn{\Delta F}).
#' @param include_species Logical. If TRUE, attempts to return a limiting-species plot as a second panel.
#' @param include_yield_share Logical. If TRUE, attempts to return an implied yield-share plot as a third panel.
#' @param title Character. Main plot title.
#' @param xlab Character. Optional x-axis label. If \code{NULL}, the function uses the species names
#'   stored in \code{attr(sens_uncon, "species")} when available.
#'
#' @return If \code{include_species = FALSE} and \code{include_yield_share = FALSE},
#' returns a single ggplot object (the main plot).
#'
#' Otherwise returns a named list with:
#' \describe{
#'   \item{\code{main}}{The main reference-point sensitivity plot (ggplot).}
#'   \item{\code{limiting_species}}{A limiting-species panel (ggplot) or \code{NULL}.}
#'   \item{\code{yield_share}}{An implied yield-share panel (ggplot) or \code{NULL}.}
#' }
#'
#' @examples
#' \dontrun{
#' sens_uncon <- spr_sensitivity(inp, prop = seq(0.2, 0.8, 0.05),
#'                              multispecies_constraint = "none")
#' sens_con   <- spr_sensitivity(inp, prop = seq(0.2, 0.8, 0.05),
#'                              multispecies_constraint = "all_species")
#'
#' # Main plot only
#' plot_recprop_sensitivity(sens_uncon, sens_con)
#'
#' # Delta plot
#' plot_recprop_sensitivity(sens_uncon, sens_con, show_delta = TRUE)
#'
#' # Return diagnostic panels (if columns are available)
#' plots <- plot_recprop_sensitivity(
#'   sens_uncon, sens_con,
#'   include_species = TRUE,
#'   include_yield_share = TRUE
#' )
#' plots$main
#' plots$limiting_species
#' plots$yield_share
#' }
#'
#' @export
plot_recprop_sensitivity <- function(
    sens_uncon,
    sens_con,
    targets = c(0.40, 0.35),
    prop_col = "prop_value",
    show_delta = FALSE,
    include_species = FALSE,
    include_yield_share = FALSE,
    title = "Sensitivity of multispecies reference points to recruitment proportions",
    xlab = NULL
) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")

  if (!is.data.frame(sens_uncon) || !is.data.frame(sens_con)) stop("sens_uncon and sens_con must be data.frames.")
  if (!prop_col %in% names(sens_uncon)) stop("prop_col not found in sens_uncon: ", prop_col)
  if (!prop_col %in% names(sens_con)) stop("prop_col not found in sens_con: ", prop_col)

  sp <- attr(sens_uncon, "species", exact = TRUE)
  if (is.null(sp) || length(sp) < 2) sp <- attr(sens_con, "species", exact = TRUE)
  if (is.null(sp) || length(sp) < 2) sp <- c("sp1", "sp2")

  if (is.null(xlab)) {
    xlab <- paste0("R0 proportion assigned to ", sp[1], " (", sp[2], " = 1 - ", sp[1], ")")
  }

  # Builds long datafrm for the different ref pt methods
  make_ref_df <- function(df, method_label, use_constrained_cols) {
    out <- lapply(targets, function(t) {
      t_int <- round(100 * t)
      col <- if (use_constrained_cols) paste0("F", t_int, "_total_constrained") else paste0("F", t_int, "_total")
      if (!col %in% names(df)) return(NULL)
      data.frame(
        prop_value = as.numeric(df[[prop_col]]),
        target = as.numeric(t),
        method = method_label,
        F = as.numeric(df[[col]]),
        stringsAsFactors = FALSE
      )
    })
    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0) return(data.frame())
    do.call(rbind, out)
  }

  df_u <- make_ref_df(sens_uncon, "Unconstrained", use_constrained_cols = FALSE)
  df_c <- make_ref_df(sens_con, "Constrained", use_constrained_cols = TRUE)
  df_long <- rbind(df_u, df_c)

  df_long <- df_long[
    is.finite(df_long$prop_value) & is.finite(df_long$target) & is.finite(df_long$F),
    ,
    drop = FALSE
  ]

  if (nrow(df_long) == 0) stop("No finite rows found to plot. Check prop_col and target columns.")

  # collapse dups
  df_long <- df_long |>
    dplyr::summarise(
      F = mean(.data$F, na.rm = TRUE),
      .by = c("prop_value", "target", "method")
    )

  # Delta mode: change in ref pt from one method to another (method - unconstrained)
  if (isTRUE(show_delta)) {
    wide <- tidyr::pivot_wider(
      df_long,
      id_cols = c("prop_value", "target"),
      names_from = "method",
      values_from = "F",
      values_fn = function(x) mean(as.numeric(x), na.rm = TRUE)
    )
    wide$Unconstrained <- as.numeric(wide$Unconstrained)
    wide$Constrained <- as.numeric(wide$Constrained)

    df_long <- rbind(
      data.frame(prop_value = wide$prop_value, target = wide$target, method = "Unconstrained", F = 0),
      data.frame(prop_value = wide$prop_value, target = wide$target, method = "Constrained",
                 F = wide$Constrained - wide$Unconstrained)
    )
    df_long <- df_long[is.finite(df_long$F), , drop = FALSE]
  }

  df_long$target <- factor(
    df_long$target,
    levels = sort(unique(df_long$target)),
    labels = paste0("SPR ", format(sort(unique(as.numeric(df_long$target))), nsmall = 2))
  )
  df_long$method <- factor(df_long$method, levels = c("Unconstrained", "Constrained"))

  ylab <- if (show_delta) "Difference in reference point F (method-unconstrained)" else "Reference point fishing mortality (F)"

  # ---- Main plot of reference point F as a function of recruitment proportion
  p_main <- ggplot2::ggplot(df_long, ggplot2::aes(x = prop_value, y = F, color = target, linetype = method)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 3) +   # 50/50 prop ref
    ggplot2::labs(x = xlab, y = ylab, color = NULL, linetype = NULL, title = title) +
    scale_color_tier4() + theme_tier4()

  if (isTRUE(show_delta)) {
    p_main <- p_main + ggplot2::geom_hline(yintercept = 0, linetype = 3) # deltaF = 0 ref
  }

  if (!isTRUE(include_species) && !isTRUE(include_yield_share)) return(p_main)

  out <- list(main = p_main)

  # ---- Optional limiting species plot ----
  # uses sens_con columns, will return NULL if missing
  if (isTRUE(include_species)) {
    lim_rows <- list()
    for (t in targets) {
      t_int <- round(100 * t)
      cols <- paste0("F", t_int, "_", sp)
      if (any(!cols %in% names(sens_con))) next

      tmp <- data.frame(
        prop_value = as.numeric(sens_con[[prop_col]]),
        target = as.numeric(t),
        F1 = as.numeric(sens_con[[cols[1]]]),
        F2 = as.numeric(sens_con[[cols[2]]]),
        stringsAsFactors = FALSE
      )
      tmp <- tmp[is.finite(tmp$prop_value) & (is.finite(tmp$F1) | is.finite(tmp$F2)), , drop = FALSE]
      tmp$limiting_species <- ifelse(tmp$F1 <= tmp$F2, sp[1], sp[2])
      lim_rows[[as.character(t)]] <- tmp[, c("prop_value","target","limiting_species")]
    }

    if (length(lim_rows) > 0) {
      lim_df <- do.call(rbind, lim_rows) |>
        dplyr::summarise(
          limiting_species = dplyr::first(.data$limiting_species),
          .by = c("prop_value", "target")
        )

      lim_df$target <- factor(
        lim_df$target,
        levels = sort(unique(lim_df$target)),
        labels = paste0("SPR ", format(sort(unique(as.numeric(lim_df$target))), nsmall = 2))
      )

      out$limiting_species <- ggplot2::ggplot(lim_df, ggplot2::aes(x = prop_value, y = 1, fill = limiting_species)) +
        ggplot2::geom_tile() +
        ggplot2::geom_vline(xintercept = 0.5, linetype = 3) +
        ggplot2::facet_wrap(~ target, ncol = 1) +
        ggplot2::scale_y_continuous(breaks = NULL) +
        ggplot2::labs(x = xlab, y = NULL, fill = "Limiting species",
                      title = "Limiting species at target (smallest species-specific F)") +
        scale_color_tier4() + theme_tier4()
    } else {
      out$limiting_species <- NULL
    }
  }

  # ---- Optional yield share plot ----
  if (isTRUE(include_yield_share)) {
    share_rows <- list()
    for (t in targets) {
      t_int <- round(100 * t)

      cols_u <- paste0("share", t_int, "_", sp)
      if (all(cols_u %in% names(sens_uncon))) {
        share_rows[[paste0("u_", t_int)]] <- rbind(
          data.frame(prop_value = as.numeric(sens_uncon[[prop_col]]), target = as.numeric(t), method = "Unconstrained",
                     species = sp[1], share = as.numeric(sens_uncon[[cols_u[1]]]), stringsAsFactors = FALSE),
          data.frame(prop_value = as.numeric(sens_uncon[[prop_col]]), target = as.numeric(t), method = "Unconstrained",
                     species = sp[2], share = as.numeric(sens_uncon[[cols_u[2]]]), stringsAsFactors = FALSE)
        )
      }

      cols_c <- paste0("share", t_int, "_", sp, "_constrained")
      if (all(cols_c %in% names(sens_con))) {
        share_rows[[paste0("c_", t_int)]] <- rbind(
          data.frame(prop_value = as.numeric(sens_con[[prop_col]]), target = as.numeric(t), method = "Constrained",
                     species = sp[1], share = as.numeric(sens_con[[cols_c[1]]]), stringsAsFactors = FALSE),
          data.frame(prop_value = as.numeric(sens_con[[prop_col]]), target = as.numeric(t), method = "Constrained",
                     species = sp[2], share = as.numeric(sens_con[[cols_c[2]]]), stringsAsFactors = FALSE)
        )
      }
    }

    if (length(share_rows) > 0) {
      sh <- do.call(rbind, share_rows)
      sh <- sh[is.finite(sh$prop_value) & is.finite(sh$share), , drop = FALSE]

      sh <- sh |>
        dplyr::summarise(
          share = mean(.data$share, na.rm = TRUE),
          .by = c("prop_value", "target", "method", "species")
        )

      sh$target <- factor(
        sh$target,
        levels = sort(unique(sh$target)),
        labels = paste0("SPR ", format(sort(unique(as.numeric(sh$target))), nsmall = 2))
      )
      sh$method <- factor(sh$method, levels = c("Unconstrained", "Constrained"))

      out$yield_share <- ggplot2::ggplot(sh |>
                                           dplyr::filter(species == sp[1]),
                                         ggplot2::aes(x = prop_value, y = share, color = target, linetype = method)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_vline(xintercept = 0.5, linetype = 3) +
        ggplot2::geom_hline(yintercept = 0.5, linetype = 3) +
        ggplot2::facet_wrap(~ species, ncol = 1) +
        ggplot2::coord_cartesian(ylim = c(0, 1)) +
        ggplot2::labs(x = xlab, y = paste0("Yield share of ", sp[1]),
                      color = NULL,
                      linetype = NULL,
                      title = paste0("Implied yield share of ", sp[1], " at reference points")) +
        scale_color_tier4() + theme_tier4()
    } else {
      out$yield_share <- NULL
    }
  }

  out
}

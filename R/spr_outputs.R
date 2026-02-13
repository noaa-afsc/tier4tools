
#' Plot SPR–F curves with reference points
#'
#' Visualize spawning potential ratio (SPR) as a function of fishing mortality (F)
#' from a \code{run_spr()} output, with report-ready annotations of SPR-based fishing
#' mortality reference points (for example, F40% and F35%).
#'
#' @details
#' \strong{What is plotted.} The function plots \code{SPR(F) = SBPR(F) /
#' SBPR(0)} across the fishing mortality grid \code{spr_out$F_vec}. Horizontal
#' reference lines are drawn at each target SPR level in
#' \code{spr_out$spr_targets} and vertical reference lines are drawn at each
#' F-based reference points.
#'
#' \strong{Single species.} If \code{spr_out} contains only one species, the function always
#' produces the single-species plot and ignores \code{which}.
#'
#' \strong{Multispecies complexes.}
#' \itemize{
#'   \item \code{which = "species"} shows one SPR(F) curve per species (colored) and annotates
#'   each species' own reference points.
#'   \item \code{which = "total"} shows the combined complex SPR(F) curve (black) with optional
#'   species overlays (colored). If constrained complex reference points exist (i.e.,
#'   \code{spr_out$F_spr_total_constrained} contains finite values), the plot also shows both
#'   unconstrained and constrained complex reference points, and highlights the limiting species
#'   at the constrained reference points.
#' }
#'
#' \strong{Multispecies constraints.} If \code{run_spr()} was called with
#' \code{multispecies_constraint = "all_species"}, constrained combined reference points may be
#' available in \code{spr_out$F_spr_total_constrained}. When present, constrained and unconstrained
#' reference points are shown together and distinguished by \code{Method} (shape and linetype).
#'
#' @param spr_out Output list from \code{run_spr()}.
#' @param which Character string. For multispecies outputs only:
#'   \describe{
#'     \item{\code{"total"}}{Plot combined SPR(F) for the complex (black), with species overlays (colored).}
#'     \item{\code{"species"}}{Plot species-specific SPR(F) curves (colored).}
#'   }
#' @param digits Integer number of digits to show in reported F values.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' spr_out1 <- run_spr(inp_single, diagnostics = FALSE)
#' plot_spr_curves(spr_out1)
#'
#' spr_out2 <- run_spr(inp_multi, diagnostics = FALSE, multispecies_constraint = "all_species")
#' plot_spr_curves(spr_out2, which = "total")
#' plot_spr_curves(spr_out2, which = "species")
#' }
#'
#' @export
plot_spr_curves <- function(
    spr_out,
    which = c("total", "species"),
    digits = 3
) {
  which <- match.arg(which)

  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")
  if (!requireNamespace("tidyr", quietly = TRUE)) stop("Please install tidyr.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install dplyr.")

  if (!is.list(spr_out) || is.null(spr_out$F_vec) || is.null(spr_out$spr_targets) || is.null(spr_out$SPR_total)) {
    stop("spr_out must be the output of run_spr().")
  }

  F <- spr_out$F_vec
  targets <- spr_out$spr_targets
  t_int <- round(100 * targets)
  fmt <- function(x) formatC(x, format = "f", digits = digits)
  idx_near <- function(Fval) which.min(abs(F - Fval))

  has_species <- !is.null(spr_out$SPR_by_species) && length(spr_out$SPR_by_species) > 0
  sp_names <- if (has_species) names(spr_out$SPR_by_species) else character(0)
  nsp <- length(sp_names)

  # extract named numeric values
  get_named_vals <- function(x, nms) {
    if (is.null(x) || !is.numeric(x)) return(setNames(rep(NA_real_, length(nms)), nms))
    out <- as.numeric(x[nms])
    names(out) <- nms
    out
  }

  # ref point marks for total SPR curve
  make_marks_total <- function(Fvals_named, method_label) {
    ok <- is.finite(Fvals_named)
    if (!any(ok)) return(data.frame())
    Fv <- as.numeric(Fvals_named[ok])
    i <- vapply(Fv, idx_near, integer(1))
    data.frame(
      target = names(Fvals_named)[ok],
      method = method_label,
      F = Fv,
      SPR = as.numeric(spr_out$SPR_total[i]),
      stringsAsFactors = FALSE
    )
  }

  # ref point marks for species SPR curves
  make_marks_species <- function(Fvals_named_by_sp) {
    # returns one row per (species, target) with point at the species curve and a segment down to x-axis
    out <- list()
    for (sp in names(Fvals_named_by_sp)) {
      Fv_named <- Fvals_named_by_sp[[sp]]
      ok <- is.finite(Fv_named)
      if (!any(ok)) next
      Fv <- as.numeric(Fv_named[ok])
      i <- vapply(Fv, idx_near, integer(1))
      out[[sp]] <- data.frame(
        species = sp,
        target = names(Fv_named)[ok],
        F = Fv,
        SPR = as.numeric(spr_out$SPR_by_species[[sp]][i]),
        stringsAsFactors = FALSE
      )
    }
    dplyr::bind_rows(out)
  }

  # ids the limiting species at constrained F
  limiting_species <- function(Fval, target) {
    if (!has_species || !is.finite(Fval)) return(NA_character_)
    i <- idx_near(Fval)
    spr_at <- vapply(sp_names, function(s) as.numeric(spr_out$SPR_by_species[[s]][i]), numeric(1))
    # choose species with smallest (SPR - target); most "binding"
    sp_names[which.min(spr_at - target)]
  }

  # extract F targets
  nm_total_un <- paste0("F", t_int, "_total")
  F_un <- get_named_vals(spr_out$F_spr_total, nm_total_un)
  names(F_un) <- paste0("F", t_int)

  nm_total_con <- paste0("F", t_int, "_total_constrained")
  F_con_raw <- get_named_vals(spr_out$F_spr_total_constrained, nm_total_con)
  names(F_con_raw) <- paste0("F", t_int)

  has_constrained <- any(is.finite(F_con_raw))

  # ---- Single species plot ----

  # for single spp analyses, always show single curve with the unconstrained
  # total reference points
  if (!has_species || nsp <= 1) {
    df <- data.frame(F = F, SPR = spr_out$SPR_total)

    marks <- make_marks_total(F_un, "Unconstrained")

    subtitle <- NULL
    ok <- is.finite(F_un)
    if (any(ok)) subtitle <- paste0(names(F_un)[ok], "=", fmt(F_un[ok]), collapse = ", ")

    p <- ggplot2::ggplot(df, ggplot2::aes(F, SPR)) +
      ggplot2::geom_line(linewidth = 1.1) +
      ggplot2::geom_hline(yintercept = targets, linetype = 2) +
      ggplot2::labs(
        x = "Fishing mortality (F)",
        y = "SPR",
        subtitle = subtitle
      )

    if (nrow(marks) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data = marks,
          ggplot2::aes(x = F, xend = F, y = SPR, yend = 0),
          inherit.aes = FALSE,
          linewidth = 0.7
        ) +
        ggplot2::geom_point(
          data = marks,
          ggplot2::aes(x = F, y = SPR),
          inherit.aes = FALSE,
          size = 2.5
        )
    }

    return(p + theme_tier4())
  }

  # ---- Multispecies: species plot ----
  if (which == "species") {
    df_sp <- data.frame(F = F, spr_out$SPR_by_species, check.names = FALSE)
    df_long <- tidyr::pivot_longer(df_sp, cols = -F, names_to = "species", values_to = "SPR")

    # species-specific F targets (from spr_out$F_spr_by_species)
    # expect names like F40_RE, F35_RE
    Ft_by_sp <- lapply(sp_names, function(sp) {
      ft <- spr_out$F_spr_by_species[[sp]]
      want <- paste0("F", t_int, "_", sp)
      v <- as.numeric(ft[want])
      names(v) <- paste0("F", t_int)
      v
    })
    names(Ft_by_sp) <- sp_names

    marks_sp <- make_marks_species(Ft_by_sp)

    # subtitle: "RE: F40=..., F35=... | BS: ..."
    parts <- character(0)
    for (sp in sp_names) {
      v <- Ft_by_sp[[sp]]
      ok <- is.finite(v)
      if (any(ok)) parts <- c(parts, paste0(sp, ": ", paste0(names(v)[ok], "=", fmt(v[ok]), collapse = ", ")))
    }
    subtitle <- if (length(parts) == 0) NULL else paste(parts, collapse = "  |  ")

    p <- ggplot2::ggplot(df_long, ggplot2::aes(F, SPR, color = species)) +
      ggplot2::geom_line(linewidth = 1.0) +
      ggplot2::geom_hline(yintercept = targets, linetype = 2) +
      ggplot2::labs(
        x = "Fishing mortality (F)",
        y = "SPR",
        color = "Species",
        subtitle = subtitle
      )

    if (nrow(marks_sp) > 0) {
      p <- p +
        ggplot2::geom_segment(
          data = marks_sp,
          ggplot2::aes(x = F, xend = F, y = SPR, yend = 0, color = species),
          inherit.aes = FALSE,
          linewidth = 0.7,
          show.legend = FALSE
        ) +
        ggplot2::geom_point(
          data = marks_sp,
          ggplot2::aes(x = F, y = SPR, color = species),
          inherit.aes = FALSE,
          size = 2.2,
          show.legend = FALSE
        )
    }

    return(p + scale_color_tier4() + theme_tier4())
  }

  # ---- Multispecies: total complex plot ----
  df_total <- data.frame(F = F, SPR = spr_out$SPR_total)
  df_sp <- data.frame(F = F, spr_out$SPR_by_species, check.names = FALSE)
  df_sp_long <- tidyr::pivot_longer(df_sp, cols = -F, names_to = "species", values_to = "SPR")

  # marks: always include unconstrained, only include constrained if it exists
  marks_un <- make_marks_total(F_un, "Unconstrained")
  marks_con <- if (has_constrained) make_marks_total(F_con_raw, "Constrained") else data.frame()
  marks <- dplyr::bind_rows(marks_un, marks_con)
  if (nrow(marks) > 0) marks$method <- factor(marks$method, levels = c("Unconstrained", "Constrained"))

  # limiting-species points at constrained F
  lim_pts <- data.frame()
  if (has_constrained) {
    ok <- is.finite(F_con_raw)
    if (any(ok)) {
      Fv <- as.numeric(F_con_raw[ok])
      targ <- targets[ok]
      lim <- mapply(limiting_species, Fv, targ, USE.NAMES = FALSE)
      ii <- vapply(Fv, idx_near, integer(1))
      lim_spr <- mapply(function(sp, i) as.numeric(spr_out$SPR_by_species[[sp]][i]), lim, ii, USE.NAMES = FALSE)

      lim_pts <- data.frame(
        method = factor("Constrained", levels = c("Unconstrained", "Constrained")),
        limiting = lim,
        F = Fv,
        SPR = as.numeric(lim_spr),
        stringsAsFactors = FALSE
      )
    }
  }

  # subtitle
  subtitle_parts <- character(0)
  ok_un <- is.finite(F_un)
  if (any(ok_un)) {
    subtitle_parts <- c(
      subtitle_parts,
      paste0("Unconstrained: ", paste0(names(F_un)[ok_un], "=", fmt(F_un[ok_un]), collapse = ", "))
    )
  }
  if (has_constrained) {
    okc <- is.finite(F_con_raw)
    if (any(okc)) {
      lim <- mapply(limiting_species, as.numeric(F_con_raw[okc]), targets[okc], USE.NAMES = FALSE)
      # "equiv complex Fspr=Fxx" uses total SPR at that F
      i <- vapply(as.numeric(F_con_raw[okc]), idx_near, integer(1))
      spr_tot_at <- as.numeric(spr_out$SPR_total[i])
      eq <- paste0("F", round(100 * spr_tot_at))

      pieces <- paste0(
        names(F_con_raw)[okc], "=", fmt(as.numeric(F_con_raw[okc])),
        ifelse(is.na(lim), "", paste0(" (", lim, " limiting), ")),
        "equiv complex Fspr=", eq
      )
      subtitle_parts <- c(subtitle_parts, paste0("Constrained: ", paste(pieces, collapse = "\n")))
    } else {
      subtitle_parts <- c(subtitle_parts, "Constrained: target not met on F grid for one or more targets")
    }
  }
  subtitle <- if (length(subtitle_parts) == 0) NULL else paste(subtitle_parts, collapse = "\n")

  # plot
  p <- ggplot2::ggplot() +
    # species curves first (colored)
    ggplot2::geom_line(
      data = df_sp_long,
      ggplot2::aes(F, SPR, color = species),
      linewidth = 0.8,
      alpha = 0.85
    ) +
    # total curve on top (black)
    ggplot2::geom_line(
      data = df_total,
      ggplot2::aes(F, SPR),
      linewidth = 1.2
    ) +
    ggplot2::geom_hline(yintercept = targets, linetype = 2) +
    ggplot2::labs(
      x = "Fishing mortality (F)",
      y = "SPR",
      color = "Species",
      subtitle = subtitle
    )

  # add marks
  if (nrow(marks) > 0) {
    if (has_constrained) {
      # both methods: map linetype + shape to method
      p <- p +
        ggplot2::geom_segment(
          data = marks,
          ggplot2::aes(x = F, xend = F, y = SPR, yend = 0, linetype = method),
          inherit.aes = FALSE,
          linewidth = 0.7
        ) +
        ggplot2::geom_point(
          data = marks,
          ggplot2::aes(x = F, y = SPR, shape = method),
          inherit.aes = FALSE,
          size = 2.4
        ) +
        ggplot2::scale_linetype_discrete(drop = FALSE, name = "Method") +
        ggplot2::scale_shape_discrete(drop = FALSE, name = "Method")
    } else {
      # unconstrained only: no Method mapping
      p <- p +
        ggplot2::geom_segment(
          data = marks,
          ggplot2::aes(x = F, xend = F, y = SPR, yend = 0),
          inherit.aes = FALSE,
          linewidth = 0.7
        ) +
        ggplot2::geom_point(
          data = marks,
          ggplot2::aes(x = F, y = SPR),
          inherit.aes = FALSE,
          size = 2.4
        )
    }
  }

  # limiting species: mapped to method, colored by limiting species
  if (nrow(lim_pts) > 0) {
    p <- p +
      ggplot2::geom_point(
        data = lim_pts,
        ggplot2::aes(x = F, y = SPR, shape = method, color = limiting),
        inherit.aes = FALSE,
        size = 2.8,
        stroke = 0.9,
        show.legend = FALSE
      )
  }

  p + scale_color_tier4() + theme_tier4()
}


#' Plot SPR diagnostics, SBPR age decompositions
#'
#' Creates report-ready diagnostic plots of age-specific decompositions of
#' spawning biomass per recruit (SBPR) under unfished conditions (F = 0) and at
#' a reference fishing mortality (F = Ftarget) chosen by \code{run_spr()}
#' (default = F40).
#'
#'
#' Returned object is a named list of ggplot objects:
#' \itemize{
#'   \item \code{plots$contrib}: SBPR contribution by age, comparing F0 vs Ftarget.
#'   \item \code{plots$removed}: SBPR removed by fishing at each age (F0 - Ftarget).
#'   \item \code{plots$removed_prop}: (not shown by default) SBPR removed by
#'   fishing at each age, normalized to the species' total unfished SBPR (F0).
#'   Useful for multispecies applications to compare fractional losses among
#'   species when absolute SBPR magnitudes differ.
#'   \item \code{plots$survivorship}: survivorship by age (per recruit), F0 vs Ftarget.
#' }
#'
#' Plus-group treatment can strongly influence SBPR, especially for long-lived
#' species. When \code{drop_plus_from_plot = TRUE} (default), the terminal age
#' is excluded from the plots to avoid dominance by the plus group, while the
#' fraction of total SBPR attributable to the terminal age is retained in
#' \code{spr_out$diagnostics$plus_share} and can be appended as a plot caption.
#'
#' \strong{Summary table.} The function also attaches a compact per-species summary table
#' (unfished SBPR, SBPR at Ftarget, implied SPR at Ftarget, and plus-group shares when available)
#' as \code{attr(plots, "summary")}.
#'
#' @param spr_out Output list from \code{run_spr()} with \code{diagnostics = TRUE}.
#' @param species Character vector of species to include. Defaults to all
#'   species in \code{spr_out}.
#' @param drop_plus_from_plot Logical. If TRUE (default), drops the terminal age from the
#'   plotted panels (only when plus-group diagnostics are available).
#' @param include_plus_caption Logical. If TRUE (default), adds the plus-group share
#'   caption to \code{plots$contrib} when available.
#' @param which_panels Character vector, any of \code{c("contrib", "removed",
#'   "removed_prop", "survivorship")}.
#'
#' @return A named list of ggplot objects. A compact summary table is attached as
#'   \code{attr(out, "summary")}.
#'
#' @examples
#' \dontrun{
#' spr_out <- run_spr(inp, diagnostics = TRUE,
#'                   multispecies_constraint = "all_species")
#'
#' # Default panels
#' plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
#' plots$contrib
#' plots$removed
#' plots$survivorship
#'
#' # Summary table for optional reporting
#' attr(plots, "summary")
#' }
#'
#' @export
plot_spr_decomp <- function(
    spr_out,
    species = NULL,
    drop_plus_from_plot = TRUE,
    include_plus_caption = TRUE,
    which_panels = c("contrib", "removed", "survivorship")
) {
  if (!is.list(spr_out) || is.null(spr_out$diagnostics) || is.null(spr_out$diagnostics$decomp)) {
    stop("spr_out must be the output of run_spr(diagnostics=TRUE).")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install ggplot2.")

  df <- spr_out$diagnostics$decomp

  # Determine available species
  all_species <- unique(df$species)
  if (is.null(species)) species <- all_species
  if (any(!species %in% all_species)) {
    stop("species must be a subset of: ", paste(all_species, collapse = ", "))
  }
  df <- df[df$species %in% species, , drop = FALSE]

  # Plotting mode
  plot_species <- unique(df$species)
  multispecies_mode <- length(plot_species) > 1

  # Ftarget = primary solid line; F0 = dashed.
  if (!all(c("F0", "Ftarget") %in% unique(df$Flevel))) {
    # allow other labels, but still coerce if these exist
    df$Flevel <- as.character(df$Flevel)
  } else {
    df$Flevel <- factor(df$Flevel, levels = c("Ftarget", "F0"))
  }

  # Plus group diagnostics
  has_plus <- !is.null(spr_out$diagnostics$plus_share) &&
    is.data.frame(spr_out$diagnostics$plus_share) &&
    nrow(spr_out$diagnostics$plus_share) > 0

  # Optionally drop plus group
  if (isTRUE(drop_plus_from_plot) && isTRUE(has_plus)) {
    last_age <- max(df$age, na.rm = TRUE)
    df_plot <- df[df$age < last_age, , drop = FALSE]
  } else {
    df_plot <- df
  }

  # Caption
  cap <- NULL
  if (isTRUE(include_plus_caption) && isTRUE(has_plus)) {
    ps <- spr_out$diagnostics$plus_share
    ps <- ps[ps$species %in% species, , drop = FALSE]
    by_sp <- split(ps, ps$species)

    parts <- vapply(names(by_sp), function(spn) {
      tmp <- by_sp[[spn]]
      p0  <- tmp$plus_share[tmp$Flevel == "F0"][1]
      pF  <- tmp$plus_share[tmp$Flevel == "Ftarget"][1]
      if (multispecies_mode) {
        sprintf("%s plus-share: F0=%s, Ftarget=%s", spn,
                sprintf("%.0f%%", 100 * p0), sprintf("%.0f%%", 100 * pF))
      } else {
        sprintf("Plus-group share: F0=%s, Ftarget=%s",
                sprintf("%.0f%%", 100 * p0), sprintf("%.0f%%", 100 * pF))
      }
    }, character(1))

    cap <- paste(parts, collapse = "   |   ")
  }

  plots <- list()

  add_defaults <- function(p, use_color_scale) {
    if (use_color_scale) p <- p + scale_color_tier4()
    p + theme_tier4()
  }

  # Single shared linetype scale for all panels that use Flevel
  add_Flevel_linetype <- function(p) {
    p + ggplot2::scale_linetype_manual(
      values = c(Ftarget = "solid", F0 = "dotted"),
      breaks = c("Ftarget", "F0"),
      labels = c("F at SPR target", "Unfished (F0)")
    )
  }

  # Summary table
  df0_full <- df[df$Flevel == "F0", , drop = FALSE]
  dfF_full <- df[df$Flevel == "Ftarget", , drop = FALSE]

  sum_by_sp <- lapply(plot_species, function(spn) {
    sb0 <- sum(df0_full$contrib[df0_full$species == spn], na.rm = TRUE)
    sbF <- sum(dfF_full$contrib[dfF_full$species == spn], na.rm = TRUE)
    spr <- if (is.finite(sb0) && sb0 > 0) sbF / sb0 else NA_real_
    data.frame(
      species = spn,
      SBPR_F0 = sb0,
      SBPR_Ftarget = sbF,
      SPR_at_Ftarget = spr,
      stringsAsFactors = FALSE
    )
  })
  summary_tbl <- do.call(rbind, sum_by_sp) # becomes attribute

  if (isTRUE(has_plus)) {
    ps <- spr_out$diagnostics$plus_share
    ps <- ps[ps$species %in% species, , drop = FALSE]
    p0 <- tapply(ps$plus_share[ps$Flevel == "F0"], ps$species[ps$Flevel == "F0"], function(z) z[1])
    pF <- tapply(ps$plus_share[ps$Flevel == "Ftarget"], ps$species[ps$Flevel == "Ftarget"], function(z) z[1])
    summary_tbl$plus_share_F0 <- as.numeric(p0[summary_tbl$species])
    summary_tbl$plus_share_Ftarget <- as.numeric(pF[summary_tbl$species])
  }

  # ---- Panel: contributions-by-age ----
  if ("contrib" %in% which_panels) {
    if (multispecies_mode) {
      p <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = age, y = contrib, color = species, linetype = Flevel)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(x = "Age", y = "SBPR contribution by age",
                      color = "Species", linetype = "Fishing level")
    } else {
      p <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = age, y = contrib, linetype = Flevel)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(x = "Age", y = "SBPR contribution by age",
                      linetype = "Fishing level")
    }

    p <- p +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
    if (!is.null(cap)) p <- p + ggplot2::labs(caption = cap)

    p <- add_Flevel_linetype(p)
    plots$contrib <- add_defaults(p, use_color_scale = multispecies_mode)
  }

  # ---- Panel: removed SBPR (absolute) ----
  if ("removed" %in% which_panels) {
    df_rem <- df_plot[df_plot$Flevel == "F0", , drop = FALSE]

    if (multispecies_mode) {
      p <- ggplot2::ggplot(df_rem, ggplot2::aes(x = age, y = removed, color = species)) +
        ggplot2::geom_line(linewidth = 1, linetype = "solid", show.legend = TRUE) +
        ggplot2::labs(x = "Age", y = "SBPR removed (F0 - Ftarget)", color = "Species")
    } else {
      p <- ggplot2::ggplot(df_rem, ggplot2::aes(x = age, y = removed)) +
        ggplot2::geom_line(linewidth = 1, linetype = "solid") +
        ggplot2::labs(x = "Age", y = "SBPR removed (F0 - Ftarget)")
    }

    p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    plots$removed <- add_defaults(p, use_color_scale = multispecies_mode)
  }

  # ---- Optional Panel: removed SBPR (proportional) ----
  if ("removed_prop" %in% which_panels) {
    df_rem <- df_plot[df_plot$Flevel == "F0", , drop = FALSE]
    denom <- tapply(df0_full$contrib, df0_full$species, sum, na.rm = TRUE)
    df_rem$removed_prop <- df_rem$removed / as.numeric(denom[df_rem$species])
    df_rem$removed_prop[!is.finite(df_rem$removed_prop)] <- NA_real_

    if (multispecies_mode) {
      p <- ggplot2::ggplot(df_rem, ggplot2::aes(x = age, y = removed_prop, color = species)) +
        ggplot2::geom_line(linewidth = 1, linetype = "solid") +
        ggplot2::labs(x = "Age", y = "SBPR removed / total SBPR(F0)", color = "Species")
    } else {
      p <- ggplot2::ggplot(df_rem, ggplot2::aes(x = age, y = removed_prop)) +
        ggplot2::geom_line(linewidth = 1, linetype = "solid") +
        ggplot2::labs(x = "Age", y = "SBPR removed / total SBPR(F0)")
    }

    p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    plots$removed_prop <- add_defaults(p, use_color_scale = multispecies_mode)
  }

  # ---- Panel: survivorship-by-age ----
  if ("survivorship" %in% which_panels) {
    if (multispecies_mode) {
      p <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = age, y = survivorship, color = species, linetype = Flevel)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(x = "Age", y = "Survivorship (per recruit)",
                      color = "Species", linetype = "Fishing level")
    } else {
      p <- ggplot2::ggplot(df_plot,
                           ggplot2::aes(x = age, y = survivorship, linetype = Flevel)
      ) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::labs(x = "Age", y = "Survivorship (per recruit)",
                      linetype = "Fishing level")
    }

    p <- p + ggplot2::theme(panel.grid = ggplot2::element_blank())
    p <- add_Flevel_linetype(p)
    plots$survivorship <- add_defaults(p, use_color_scale = multispecies_mode)
  }

  attr(plots, "summary") <- summary_tbl
  return(plots)
}

#' Arrange SPR decomposition plots into report-ready grid
#'
#' Combined list output from \code{plot_spr_decomp()} into a single multipanel
#' figure using \pkg{patchwork}. Because \pkg{patchwork} is in Suggests, this
#' function will only run when it is installed.
#'
#' @details
#' Typical usage:
#' \preformatted{
#'   decomp_plots <- plot_spr_decomp(spr_out, drop_plus_from_plot = TRUE)
#'   grid_spr_decomp(decomp_plots, order = c("contrib","removed","survivorship"), ncol = 1)
#' }
#'
#' @param plots A named list of ggplot objects returned by \code{plot_spr_decomp()}.
#' @param order Character vector giving the plot names (in \code{plots}) in the order to arrange.
#'   Defaults to \code{c("contrib","removed","removed_prop","survivorship")} when present,
#'   otherwise \code{names(plots)}.
#' @param ncol Number of columns in the grid.
#' @param nrow Optional number of rows. If \code{NULL}, patchwork chooses automatically.
#' @param guides Passed to \code{patchwork::plot_layout(guides = ...)}. Default \code{"collect"} collects legends when possible.
#' @param title Optional figure title (character). Added with \code{patchwork::plot_annotation()}.
#' @param caption Optional figure caption (character). Added with \code{patchwork::plot_annotation()}.
#'
#' @return A patchwork object (a ggplot-like object) that can be printed or saved.
#' @export
grid_spr_decomp <- function(
    plots,
    order = NULL,
    ncol = 1,
    nrow = NULL,
    guides = c("collect", "keep"),
    title = NULL,
    caption = NULL
) {
  guides <- match.arg(guides)

  if (!is.list(plots) || length(plots) == 0) {
    stop("grid_spr_decomp(): 'plots' must be a non-empty list returned by plot_spr_decomp().")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("grid_spr_decomp(): Please install 'patchwork' to arrange plots (it is in Suggests).")
  }

  default_order <- c("contrib", "removed", "removed_prop", "survivorship")
  if (is.null(order)) {
    order <- default_order[default_order %in% names(plots)]
    if (length(order) == 0) order <- names(plots)
  }

  missing <- setdiff(order, names(plots))
  if (length(missing) > 0) {
    stop("grid_spr_decomp(): requested plot(s) not found: ", paste(missing, collapse = ", "))
  }

  plot_list <- plots[order]

  pw <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol, nrow = nrow) +
    patchwork::plot_layout(guides = guides)

  if (!is.null(title) || !is.null(caption)) {
    pw <- pw + patchwork::plot_annotation(title = title, caption = caption)
  }

  pw
}

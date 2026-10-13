#!/usr/bin/env Rscript
# One panel per PDF. Run: Rscript plot_three_metrics.R input.tsv output_dir [--overwrite]
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
  library(jsonlite)
})
setDTthreads(6)
args <- commandArgs(trailingOnly = TRUE)
overwrite <- "--overwrite" %in% args
args <- args[args != "--overwrite"]
if (length(args) != 2L) {
  stop("Usage: Rscript plot_three_metrics.R three_metrics.tsv output_dir [--overwrite]", call. = FALSE)
}
input <- normalizePath(args[1], mustWork = TRUE)
out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
manifest_path <- file.path(dirname(input), "analysis_summary.json")
stopifnot(file.exists(manifest_path))
manifest <- fromJSON(manifest_path)
if (is.null(manifest$contigs)) {
  stop("This example requires paired genome-mode results with reference contig lengths.", call. = FALSE)
}
control <- manifest$sample_labels$control
test <- manifest$sample_labels$test
if (is.null(control)) control <- "control"
if (is.null(test)) test <- "test"
N_cpg <- manifest$config$k
names_in <- c(paste0(control, "_beta"), paste0(test, "_beta"), "delta_beta",
              paste0(control, "_entropy"), paste0(test, "_entropy"), "delta_entropy", "epiallele_distance")
names_work <- c("beta_control", "beta_test", "delta_beta", "entropy_control", "entropy_test", "delta_entropy", "distance")
message("Reading and validating ", input)
dt <- fread(input, select = c("chrom", "start", "end", names_in), nThread = 6)
setnames(dt, names_in, names_work)
n_units <- nrow(dt)
stopifnot(n_units > 0L, all(dt$start >= 0), all(dt$end > dt$start), !anyNA(dt$chrom))
stopifnot(all(vapply(dt[, ..names_work], function(x) all(is.finite(x)), logical(1))))
for (name in names_work) {
  lower <- if (startsWith(name, "delta")) -1 else 0
  stopifnot(all(dt[[name]] >= lower - 1e-10 & dt[[name]] <= 1 + 1e-10))
}
stopifnot(max(abs(dt$delta_beta - (dt$beta_test - dt$beta_control))) < 1e-10,
          max(abs(dt$delta_entropy - (dt$entropy_test - dt$entropy_control))) < 1e-10)
expected <- manifest$outputs[[basename(input)]]$rows
if (!is.null(expected)) stopifnot(n_units == expected)
contigs <- as.data.table(manifest$contigs)[, .(chrom, length = as.numeric(length))]
stopifnot(all(unique(dt$chrom) %in% contigs$chrom))
stopifnot(all(dt$end <= contigs$length[match(dt$chrom, contigs$chrom)]))
chrom_order <- contigs$chrom
subject <- basename(dirname(input))
footer <- sprintf("%s | %s to %s | N = %s | %s CpG units\nComplete metrics only; no minimum-depth filter applied. Exploratory analysis of one sample pair.",
                  subject, control, test, N_cpg, comma(n_units))
colors <- setNames(c("#24788C", "#CF683E"), c(control, test))
theme_set(theme_bw(base_size = 11, base_family = "DejaVu Sans") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "#E9EDF0", linewidth = 0.25),
        panel.border = element_rect(color = "#61707B", linewidth = 0.45),
        plot.title = element_text(size = 15, face = "bold", color = "#213743", margin = margin(b = 8)),
        plot.subtitle = element_text(size = 10, color = "#53636E", margin = margin(b = 12)),
        plot.caption = element_text(size = 8, color = "#596772", hjust = 0, margin = margin(t = 14)),
        axis.title = element_text(size = 11), axis.text = element_text(color = "#33444F", size = 9),
        legend.position = "top", legend.title = element_text(size = 9),
        legend.text = element_text(size = 9), plot.margin = margin(16, 18, 12, 14)))
files <- character()
descriptions <- character()
save_panel <- function(p, name, title, subtitle, width = 7.8, height = 6.0) {
  path <- file.path(out, paste0(name, ".pdf"))
  if (file.exists(path) && !overwrite) stop("Output exists; use --overwrite explicitly: ", path)
  p <- p + labs(title = title, subtitle = subtitle, caption = footer)
  ggsave(path, p, device = grDevices::cairo_pdf, width = width, height = height,
         units = "in", bg = "white", family = "DejaVu Sans")
  files <<- c(files, basename(path))
  descriptions <<- c(descriptions, paste(basename(path), title, subtitle, sep = " | "))
  message("Saved ", basename(path))
}
hist_data <- function(value, bounds = c(0, 1), nbins = 50L, label = "") {
  index <- pmax(0L, pmin(nbins - 1L, floor((value - bounds[1]) / diff(bounds) * nbins)))
  counts <- tabulate(index + 1L, nbins = nbins)
  stopifnot(sum(counts) == length(value))
  data.table(x = bounds[1] + (seq_len(nbins) - 0.5) * diff(bounds) / nbins,
             percent = 100 * counts / length(value), sample = label)
}
hist_panel <- function(columns, labels, bounds, xlabel, name, title) {
  values <- rbindlist(Map(function(col, label) hist_data(dt[[col]], bounds, label = label), columns, labels))
  if (length(columns) == 2) {
    p <- ggplot(values, aes(x, percent, color = sample)) + geom_line(linewidth = 0.9) +
      geom_point(size = 1.2) + scale_color_manual(values = colors, name = NULL)
  } else {
    p <- ggplot(values, aes(x, percent)) + geom_col(width = diff(bounds) / 50, fill = "#3D7C8B", color = NA)
  }
  p <- p + scale_x_continuous(limits = bounds, expand = expansion(mult = 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.06))) + labs(x = xlabel, y = "CpG units per bin (%)")
  save_panel(p, name, title, sprintf("All units; 50 equal-width bins (width %.2f)", diff(bounds) / 50))
}
hist_panel(c("beta_control", "beta_test"), c(control, test), c(0, 1), "Methylation level (β)",
           "01_methylation_distributions", "Methylation-level distributions")
hist_panel(c("entropy_control", "entropy_test"), c(control, test), c(0, 1), "Normalized pattern entropy (H)",
           "02_entropy_distributions", "Pattern-entropy distributions")
hist_panel("distance", "", c(0, 1), "Epiallele distance", "03_epiallele_distance_distribution", "Epiallele-distance distribution")
hist_panel("delta_beta", "", c(-1, 1), sprintf("Δβ (%s minus %s)", test, control),
           "04_delta_methylation_distribution", "Changes in methylation level")
hist_panel("delta_entropy", "", c(-1, 1), sprintf("ΔH (%s minus %s)", test, control),
           "05_delta_entropy_distribution", "Changes in pattern entropy")

bin2d <- function(x, y, bounds = c(0, 1), nbins = 80L) {
  xb <- pmax(0L, pmin(nbins - 1L, floor((x - bounds[1]) / diff(bounds) * nbins)))
  yb <- pmax(0L, pmin(nbins - 1L, floor((y - bounds[1]) / diff(bounds) * nbins)))
  grouped <- data.table(xb, yb, distance = dt$distance)[,
    .(count = .N, mean_distance = mean(distance)), by = .(xb, yb)]
  stopifnot(sum(grouped$count) == n_units)
  grouped[, `:=`(x = bounds[1] + (xb + 0.5) * diff(bounds) / nbins,
                y = bounds[1] + (yb + 0.5) * diff(bounds) / nbins)]
  grouped
}
density_panel <- function(grid, bounds, xlabel, ylabel, diagonal = FALSE) {
  p <- ggplot(grid, aes(x, y, fill = count)) +
    geom_tile(width = diff(bounds) / 80, height = diff(bounds) / 80) +
    scale_fill_viridis_c(option = "C", trans = "log10", labels = label_number(big.mark = ","), name = "CpG units") +
    scale_x_continuous(limits = bounds, expand = expansion(mult = 0)) +
    scale_y_continuous(limits = bounds, expand = expansion(mult = 0)) +
    coord_fixed() + labs(x = xlabel, y = ylabel) + theme(legend.position = "right")
  if (diagonal) p <- p + geom_abline(slope = 1, intercept = 0, color = "#283C48", linetype = "dashed", linewidth = 0.4)
  p
}
for (metric in c("beta", "entropy")) {
  x <- dt[[paste0(metric, "_control")]]; y <- dt[[paste0(metric, "_test")]]
  label <- if (metric == "beta") "Methylation level (β)" else "Pattern entropy (H)"
  id <- if (metric == "beta") "06_paired_methylation_density" else "07_paired_entropy_density"
  p <- density_panel(bin2d(x, y), c(0, 1), paste(control, label), paste(test, label), TRUE)
  save_panel(p, id, paste(control, "versus", test, if (metric == "beta") "methylation" else "entropy"),
             sprintf("All units; 80 × 80 bins; logarithmic counts; Pearson r = %s",
                     if (length(x) > 1L && sd(x) > 0 && sd(y) > 0) sprintf("%.3f", cor(x, y)) else "undefined (constant values)"))
}
joint <- bin2d(dt$delta_beta, dt$delta_entropy, c(-1, 1))
zero_lines <- list(geom_hline(yintercept = 0, linetype = "dashed", color = "#37454E", linewidth = 0.3),
                   geom_vline(xintercept = 0, linetype = "dashed", color = "#37454E", linewidth = 0.3))
p <- density_panel(joint, c(-1, 1), "Δβ (test minus control)", "ΔH (test minus control)") + zero_lines
save_panel(p, "08_joint_changes_density", "Methylation and entropy changes", "All units; 80 × 80 bins; colour shows logarithmic unit counts")
mean_grid <- joint[count >= 10]
p <- ggplot(mean_grid, aes(x, y, fill = mean_distance)) + geom_tile(width = 2 / 80, height = 2 / 80) +
  scale_fill_viridis_c(limits = c(0, 1), option = "D", name = "Mean epiallele\ndistance") +
  scale_x_continuous(limits = c(-1, 1), expand = expansion(mult = 0)) +
  scale_y_continuous(limits = c(-1, 1), expand = expansion(mult = 0)) + coord_fixed() + zero_lines +
  labs(x = "Δβ (test minus control)", y = "ΔH (test minus control)") + theme(legend.position = "right")
save_panel(p, "09_joint_changes_mean_distance", "Epiallele distance across the change map",
           sprintf("80 × 80 bins; means shown for bins with ≥10 units (%.2f%% of units)", 100 * sum(mean_grid$count) / n_units))
for (metric in c("beta", "entropy")) {
  xlabel <- if (metric == "beta") "Absolute methylation change (|Δβ|)" else "Absolute entropy change (|ΔH|)"
  id <- if (metric == "beta") "10_distance_vs_methylation_change" else "11_distance_vs_entropy_change"
  p <- density_panel(bin2d(abs(dt[[paste0("delta_", metric)]]), dt$distance), c(0, 1), xlabel, "Epiallele distance")
  save_panel(p, id, paste("Epiallele distance versus", if (metric == "beta") "methylation change" else "entropy change"),
             "All units; 80 × 80 bins; colour shows logarithmic unit counts")
}

chrom <- dt[, .(units = .N, beta_control = mean(beta_control), beta_test = mean(beta_test),
                entropy_control = mean(entropy_control), entropy_test = mean(entropy_test),
                median_distance = median(distance), q25 = quantile(distance, 0.25), q75 = quantile(distance, 0.75)), by = chrom]
chrom[, chromosome := factor(chrom, levels = rev(chrom_order))]
for (metric in c("beta", "entropy")) {
  long <- melt(chrom, id.vars = "chromosome", measure.vars = paste0(metric, c("_control", "_test")),
               variable.name = "role", value.name = "mean_value")
  long[, sample := ifelse(grepl("_control$", role), control, test)]
  p <- ggplot(long, aes(mean_value, chromosome, color = sample)) +
    geom_point(position = position_dodge(width = 0.4), size = 2) +
    scale_color_manual(values = colors, name = NULL) + scale_x_continuous(limits = c(0, 1)) +
    labs(x = if (metric == "beta") "Mean methylation level (β)" else "Mean normalized entropy (H)", y = "Reference contig")
  save_panel(p, if (metric == "beta") "12_chromosome_methylation" else "13_chromosome_entropy",
             if (metric == "beta") "Methylation by chromosome" else "Pattern entropy by chromosome",
             "Arithmetic mean across plotted CpG units; unit-weighted, not read-weighted", width = 7.8, height = 8.4)
}
p <- ggplot(chrom, aes(median_distance, chromosome)) +
  geom_segment(aes(x = q25, xend = q75, yend = chromosome), color = "#9BBCC0", linewidth = 2) +
  geom_point(color = "#285F70", size = 2.2) + scale_x_continuous(limits = c(0, 1)) +
  labs(x = "Epiallele distance", y = "Reference contig")
save_panel(p, "14_chromosome_epiallele_distance", "Epiallele distance by chromosome",
           "Points: median; segments: interquartile range across plotted units (not confidence intervals)", width = 8.5, height = 8.4)
p <- ggplot(chrom, aes(units, chromosome)) + geom_col(fill = "#468896", width = 0.65) +
  scale_x_continuous(labels = label_number(scale = 1e-6, suffix = "M"), expand = expansion(mult = c(0, 0.05))) +
  labs(x = "CpG units with complete metrics (millions)", y = "Reference contig")
save_panel(p, "15_chromosome_unit_counts", "Plotted units by chromosome",
           "Nonmissing metric rows; these counts describe data availability, not read depth", width = 7.8, height = 8.4)

message("Aggregating genomic tracks in fixed 1-Mb windows")
window_bp <- 1000000
windows <- dt[, .(units = .N, delta_beta = mean(delta_beta), delta_entropy = mean(delta_entropy),
                  distance = mean(distance)), by = .(chrom, window = floor(start / window_bp))]
contigs[, offset := c(0, head(cumsum(length), -1))]
contigs[, midpoint := (offset + length / 2) / 1e6]
windows <- merge(windows, contigs, by = "chrom", sort = FALSE)
windows[, position_mb := (offset + pmin((window + 0.5) * window_bp, length)) / 1e6]
setorder(windows, offset, window)
for (metric in c("delta_beta", "delta_entropy", "distance", "units")) {
  values <- copy(windows)
  if (metric != "units") values[units < 100, (metric) := NA_real_]
  spec <- switch(metric,
    delta_beta = c("16_genome_methylation_change", "Methylation change across the genome", "Mean Δβ (test minus control)"),
    delta_entropy = c("17_genome_entropy_change", "Entropy change across the genome", "Mean ΔH (test minus control)"),
    distance = c("18_genome_epiallele_distance", "Epiallele distance across the genome", "Mean epiallele distance"),
    units = c("19_genome_unit_counts", "Plotted-unit density across the genome", "CpG units per 1-Mb window"))
  p <- ggplot(values, aes(position_mb, .data[[metric]])) +
    geom_vline(xintercept = contigs$offset[-1] / 1e6, color = "#CCD4D9", linewidth = 0.25) +
    geom_point(color = if (metric == "distance") "#2F796C" else "#376F89", size = 0.8, alpha = 0.85, na.rm = TRUE) +
    scale_x_continuous(breaks = contigs$midpoint, labels = contigs$chrom, expand = expansion(mult = 0.008)) +
    labs(x = "Reference contigs concatenated in BAM order (positions in Mb)", y = spec[3]) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), panel.grid.major.x = element_blank())
  if (metric %in% c("delta_beta", "delta_entropy")) p <- p + geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) + scale_y_continuous(limits = c(-1, 1))
  if (metric == "distance") p <- p + scale_y_continuous(limits = c(0, 1))
  if (metric == "units") p <- p + scale_y_continuous(labels = label_number(big.mark = ","))
  save_panel(p, spec[1], spec[2], if (metric == "units")
             "All nonempty windows; units assigned by start coordinate; empty windows have no point" else
             "1-Mb windows with ≥100 plotted units; arithmetic means; omitted windows have no point",
             width = 12.5, height = 5.6)
}
notes <- c("MethylDiver exploratory plots: one panel per PDF", paste("Input:", input),
           paste("Input MD5:", unname(tools::md5sum(input))), paste("Rows:", comma(n_units)),
           paste("Samples:", control, "and", test), paste("CpGs per unit:", N_cpg),
           "All TSV rows are used. No random subsampling is performed.",
           "The TSV excludes missing metrics but has no read-depth columns; no minimum-depth filter can be applied from it alone.",
           paste("The run summary reports", manifest$n_units_passing_min_depth, "units passing its depth threshold of", manifest$config$min_depth,
                 "but the TSV does not identify those units."),
           "Overlapping units are not independent biological replicates. No significance tests or confidence intervals are plotted.",
           "Entropy and pattern distances are sensitive to complete-read depth. These figures are descriptive, not depth-controlled biological conclusions.",
           "Chromosome labels retain reference accessions exactly; no chromosome-name mapping is assumed.",
           "Histograms: 50 bins, normalized to percentage of all units. Two-dimensional density: 80 x 80 bins, log10 counts.",
           "Mean-distance change map: bins with at least 10 units; this is a display-bin count threshold, not a read-depth filter.",
           "Chromosome summaries: arithmetic sample means; median distance with interquartile ranges, not confidence intervals.",
           "Genomic tracks: fixed 1-Mb windows, assigned by unit start; metric means require at least 100 plotted units/window.",
           "Input validation: all values finite and in range; deltas verified; row count matches manifest; coordinates within reference lengths.",
           "", "Metric summaries (mean / median):",
           vapply(names_work, function(name) sprintf("%s: %.6f / %.6f", name, mean(dt[[name]]), median(dt[[name]])), character(1)),
           "", "Figure index:", descriptions, "", "R session:", capture.output(sessionInfo()))
writeLines(notes, file.path(out, "plotting_notes.txt"))
stopifnot(length(files) == 19L, all(file.info(file.path(out, files))$size > 1000))
message("Finished: ", length(files), " single-panel PDFs in ", normalizePath(out))

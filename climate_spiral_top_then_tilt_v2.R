# Climate spiral from combined_ipcc2000_nasa_monthly.csv
# Phase 1: top view, cumulative by year
# Phase 2: after final year, tilt toward side view
# Improvements:
#   1) smooth yearly loops with periodic interpolation
#   2) tilt direction puts future years toward the upper side of the image
#   3) z-scale reduced to half of the previous script
#   4) reconstructed-era monthly values can be jittered around each year's mean

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# -----------------------------
# Config
# -----------------------------
csv_path <- "./ICPP_/combined_ipcc2000_nasa_monthly.csv"
out_dir  <- "spiral_top_then_tilt_v2_frames"
mp4_path <- "climate_spiral_top_then_tilt_v2.mp4"

width_px  <- 1280
height_px <- 1280
res_dpi   <- 144
fps       <- 20

bg_col        <- "black"
ring_col      <- "#e6d200"
zero_col      <- "#00ff38"
text_col      <- ring_col
center_col    <- "#f26b8a"
center_col_tv <- "#e8f3ff"

base_radius   <- 1.75
radial_scale  <- 0.62    # radius = base_radius + radial_scale * anomaly
z_scale       <- 0.07   # half of the previous 0.028
line_width    <- 2.2
line_alpha    <- 0.42
inner_hole_r  <- 0.92

# smoothness of yearly loops
points_per_year <- 240

# reconstructed-era monthly jitter
apply_recon_jitter    <- TRUE
recon_jitter_seed     <- 20260420
recon_sigma_scale     <- 1.00   # multiply estimated monthly sigma
recon_min_sigma       <- 0.03   # avoid perfectly circular loops
recon_max_sigma       <- 0.18   # avoid excessive wobble

# top view year build
step_years    <- 1
# tilt phase
n_tilt_frames <- 10
max_tilt_deg  <- 90

# -----------------------------
# Helpers
# -----------------------------
stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("CSV file not found: ", path)
}

pick_col <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

month_names_en <- c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
month_theta <- pi/2 - 2*pi*((1:12)-1)/12

circle_xy <- function(r, n = 500) {
  th <- seq(0, 2*pi, length.out = n)
  data.frame(x = r*cos(pi/2 - th), y = r*sin(pi/2 - th))
}

col_from_anom <- function(v) {
  breaks <- c(-1.2, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4)
  cols   <- c("#4e4effaf", "#f2f5d1", "#e7f6bc", "#ebde85", "#ffb37d", "#e25f75d7")
  idx <- findInterval(v, breaks, all.inside = TRUE)
  grDevices::adjustcolor(cols[idx], alpha.f = line_alpha)
}

project_xyz_old <- function(x, y, z, tilt_deg = 0, yaw_deg = 0, perspective = 0.0015) {
  ay <- yaw_deg * pi/180
  ax <- tilt_deg * pi/180

  x1 <- x*cos(ay) - y*sin(ay)
  y1 <- x*sin(ay) + y*cos(ay)
  z1 <- z

  x2 <- x1
  y2 <- y1*cos(ax) - z1*sin(ax)
  z2 <- y1*sin(ax) + z1*cos(ax)

  s <- 1 / (1 + perspective * (z2 - min(z2, na.rm = TRUE)))
  list(X = x2 * s, Y = y2 * s, Z = z2)
}

project_xyz <- function(
  x, y, z,
  tilt_deg = 0,
  yaw_deg = 0,
  perspective = 0.0015,
  z_flatten = 0.10,
  y_flatten = 0.82
) {
  ay <- yaw_deg * pi/180
  ax <- tilt_deg * pi/180

  x1 <- x*cos(ay) - y*sin(ay)
  y1 <- x*sin(ay) + y*cos(ay)

  # Z方向を強く圧縮
  z1 <- z * z_flatten

  x2 <- x1
  y2 <- y1*cos(ax) - z1*sin(ax)
  z2 <- y1*sin(ax) + z1*cos(ax)

  s <- 1 / (1 + perspective * (z2 - min(z2, na.rm = TRUE)))

  list(
    X = x2 * s,
    Y = y2 * s * y_flatten,
    Z = z2
  )
}

open_device <- function(filename) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(filename, width = width_px, height = height_px, units = "px", res = res_dpi, background = bg_col)
  } else {
    png(filename, width = width_px, height = height_px, res = res_dpi, bg = bg_col)
  }
}

estimate_observed_monthly_sigma <- function(df_obs) {
  sigma_by_year <- df_obs %>%
    group_by(year) %>%
    summarize(sigma = stats::sd(anomaly - mean(anomaly, na.rm = TRUE), na.rm = TRUE), .groups = "drop")

  sigma <- median(sigma_by_year$sigma, na.rm = TRUE)
  if (!is.finite(sigma)) sigma <- 0.08
  sigma
}

make_recon_monthly_jitter <- function(one_year_df, sigma, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  mu <- mean(one_year_df$anomaly, na.rm = TRUE)
  n  <- nrow(one_year_df)
  vals <- rnorm(n, mean = mu, sd = sigma)

  # preserve each year's annual mean exactly
  vals <- vals - mean(vals, na.rm = TRUE) + mu
  one_year_df$anomaly <- vals
  one_year_df
}

smooth_year_loop <- function(dd, points_per_year = 240) {
  dd <- dd %>% arrange(month)

  # month 13 is the repeated January closure point
  theta_base <- dd$theta
  r_base     <- dd$r

  # periodic interpolation over [0,1]
  t_base <- seq(0, 1, length.out = length(theta_base))
  t_new  <- seq(0, 1, length.out = points_per_year)

  # use periodic interpolation on x/y to avoid angle wrap issues
  x_base <- r_base * cos(theta_base)
  y_base <- r_base * sin(theta_base)

  fx <- splinefun(t_base, x_base, method = "periodic")
  fy <- splinefun(t_base, y_base, method = "periodic")

  x_new <- fx(t_new)
  y_new <- fy(t_new)

  data.frame(
    year = dd$year[1],
    x = x_new,
    y = y_new,
    z = dd$z[1],
    anomaly_mean = mean(dd$anomaly[dd$month <= 12], na.rm = TRUE)
  )
}

# -----------------------------
# Read data
# -----------------------------
stop_if_missing(csv_path)
raw <- read_csv(csv_path, show_col_types = FALSE)

nms <- names(raw)
year_col    <- pick_col(nms, c("year"))
month_col   <- pick_col(nms, c("month"))
anom_col    <- pick_col(nms, c("anomaly_degC", "anomaly", "annual_anomaly"))
source_col  <- pick_col(nms, c("source"))

if (any(is.na(c(year_col, month_col, anom_col)))) {
  stop("Required columns not found. Need year, month, anomaly_degC/anomaly")
}

df <- raw %>%
  transmute(
    year = as.integer(.data[[year_col]]),
    month = as.integer(.data[[month_col]]),
    anomaly = as.numeric(.data[[anom_col]]),
    source = if (!is.na(source_col)) as.character(.data[[source_col]]) else ""
  ) %>%
  filter(!is.na(year), !is.na(month), !is.na(anomaly)) %>%
  filter(month >= 1, month <= 12) %>%
  arrange(year, month)

# complete months only
keep_years <- df %>% count(year) %>% filter(n >= 12) %>% pull(year)
df <- df %>% filter(year %in% keep_years)

# identify observed vs reconstructed era
if ("source" %in% names(df) && any(nzchar(df$source))) {
  src_low <- tolower(df$source)
  observed_mask <- grepl("nasa|gistemp|observ", src_low)
} else {
  observed_mask <- df$year >= 1880
}
observed_start_year <- min(df$year[observed_mask], na.rm = TRUE)
if (!is.finite(observed_start_year)) observed_start_year <- 1880L

# jitter reconstructed-era monthly values if they are flat-by-year
if (apply_recon_jitter) {
  sigma_obs <- estimate_observed_monthly_sigma(df %>% filter(year >= observed_start_year))
  sigma_use <- max(recon_min_sigma, min(recon_max_sigma, sigma_obs * recon_sigma_scale))

  set.seed(recon_jitter_seed)
  recon_years <- sort(unique(df$year[df$year < observed_start_year]))

  if (length(recon_years) > 0) {
    recon_list <- lapply(recon_years, function(yy) {
      one <- df %>% filter(year == yy)
      # only jitter years that are effectively flat within the year
      if (stats::sd(one$anomaly, na.rm = TRUE) < 1e-10) {
        make_recon_monthly_jitter(one, sigma = sigma_use)
      } else {
        one
      }
    })

    df <- bind_rows(
      bind_rows(recon_list),
      df %>% filter(year >= observed_start_year)
    ) %>% arrange(year, month)
  }
}

# close each year loop by repeating Jan at month 13
plot_df <- bind_rows(
  df,
  df %>% filter(month == 1) %>% mutate(month = 13)
) %>%
  mutate(
    theta = pi/2 - 2*pi*(month - 1)/12,
    r = base_radius + radial_scale * anomaly
  )

years <- sort(unique(plot_df$year))

y0 <- min(years)
plot_df <- plot_df %>% mutate(z = (year - y0) * z_scale)

year_summary <- df %>% group_by(year) %>% summarize(anom_mean = mean(anomaly, na.rm = TRUE), .groups = "drop")
year_cols <- setNames(vapply(year_summary$anom_mean, col_from_anom, character(1)), year_summary$year)

# smooth coordinates year by year
split_years_raw <- split(plot_df, plot_df$year)
split_years <- lapply(split_years_raw, smooth_year_loop, points_per_year = points_per_year)

all_year_df <- bind_rows(split_years)

# -----------------------------
# Layout extents
# -----------------------------
r_outer <- base_radius + radial_scale * 1.15
r_zero  <- base_radius
r_inner <- base_radius + radial_scale * (-1.0)
label_r <- r_outer + 0.38
xy_lim_top <- c(-(label_r + 0.55), label_r + 0.55)

all_x <- all_year_df$x
all_y <- all_year_df$y
all_z <- all_year_df$z
cand_tilt <- seq(0, max_tilt_deg, length.out = n_tilt_frames)
proj_ranges <- lapply(cand_tilt, function(td) {
  # negative tilt so future/high-z rises toward the upper side of the image
  pr <- project_xyz(all_x, all_y, all_z, tilt_deg = -td, yaw_deg = 0)
  c(range(pr$X, na.rm = TRUE), range(pr$Y, na.rm = TRUE))
})
proj_ranges <- do.call(rbind, proj_ranges)
xy_lim_tilt_x <- range(proj_ranges[,1:2])
xy_lim_tilt_y <- range(proj_ranges[,3:4])
pad_x <- diff(xy_lim_tilt_x) * 0.16
pad_y <- diff(xy_lim_tilt_y) * 0.16
xy_lim_tilt_x <- xy_lim_tilt_x + c(-pad_x, pad_x)
xy_lim_tilt_y <- xy_lim_tilt_y + c(-pad_y, pad_y)

# -----------------------------
# Static overlays for top view
# -----------------------------
draw_top_overlay <- function(current_year) {
  for (rr in c(r_outer, r_zero, r_inner, inner_hole_r)) {
    cc <- circle_xy(rr)
    lines(cc$x, cc$y, col = if (abs(rr - r_zero) < 1e-9) zero_col else ring_col, lwd = if (abs(rr - r_zero) < 1e-9) 2.5 else 2.0)
  }

  tx <- label_r * cos(month_theta)
  ty <- label_r * sin(month_theta)
  text(tx, ty, labels = month_names_en, col = text_col, cex = 1.2, font = 2)

  text(0, r_outer, "+1°C", col = text_col, cex = 1.2, pos = 3, font = 2)
  text(0, r_zero,  "0°",   col = zero_col, cex = 1.2, pos = 3, font = 2)
  text(0, r_inner, "-1°C", col = text_col, cex = 1.2, pos = 3, font = 2)

  text(0, 0, labels = sprintf("A.D.%03d", current_year), col = center_col_tv, cex = 2.2, font = 2)
}

# -----------------------------
# Draw phases
# -----------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
frame_id <- 1L

# Phase 1: top view, cumulative by year
build_years <- years[seq(1, length(years), by = step_years)]
observed_start_year <- 1880

build_years <- c(
  years[years < 50],
  years[years >= 50 & years < observed_start_year & years %% 100 == 0],
  years[years >= observed_start_year]
)

build_years <- sort(unique(build_years))

if (tail(build_years, 1) != tail(years, 1)) build_years <- c(build_years, tail(years, 1))

for (yy in build_years) {
  fn <- file.path(out_dir, sprintf("frame-%06d.png", frame_id))
  open_device(fn)
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = xy_lim_top, ylim = xy_lim_top, asp = 1)

  for (ycur in years[years <= yy]) {
    dd <- split_years[[as.character(ycur)]]
    lines(dd$x, dd$y, col = year_cols[as.character(ycur)], lwd = line_width)
  }

  draw_top_overlay(yy)
  dev.off()
  frame_id <- frame_id + 1L
}

fps <- 20
hold_frames <- fps * 2

last_frame_idx <- frame_id - 1
last_frame_file <- file.path(out_dir, sprintf("frame-%06d.png", last_frame_idx))

for (i in seq_len(hold_frames)) {
  new_idx <- frame_id
  new_file <- file.path(out_dir, sprintf("frame-%06d.png", new_idx))
  file.copy(last_frame_file, new_file, overwrite = TRUE)
  frame_id <- frame_id + 1
}

# Phase 2: tilt after final year is fully drawn
last_year <- max(years)
ordered_years <- years
for (i in seq_len(n_tilt_frames)) {
  tilt <- max_tilt_deg * (i - 1) / max(1, n_tilt_frames - 1)
  fn <- file.path(out_dir, sprintf("frame-%06d.png", frame_id))
  open_device(fn)
  par(mar = c(0,0,0,0), xaxs = "i", yaxs = "i")
  plot.new()
  plot.window(xlim = xy_lim_tilt_x, ylim = xy_lim_tilt_y, asp = 1)

  for (ycur in ordered_years) {
    dd <- split_years[[as.character(ycur)]]
    pr <- project_xyz(dd$x, dd$y, dd$z, tilt_deg = -tilt, yaw_deg = 0)
    lines(pr$X, pr$Y, col = year_cols[as.character(ycur)], lwd = line_width)
  }

  #text(mean(xy_lim_tilt_x), quantile(xy_lim_tilt_y, 0.67), labels = sprintf("%d", last_year), col = center_col, cex = 2.0, font = 2)
  dev.off()
  frame_id <- frame_id + 1L
}

fps <- 20
hold_frames <- fps * 3

last_frame_idx <- frame_id - 1
last_frame_file <- file.path(out_dir, sprintf("frame-%06d.png", last_frame_idx))

for (i in seq_len(hold_frames)) {
  new_idx <- frame_id
  new_file <- file.path(out_dir, sprintf("frame-%06d.png", new_idx))
  file.copy(last_frame_file, new_file, overwrite = TRUE)
  frame_id <- frame_id + 1
}
message("Frames written to: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))

# optional mp4
ffmpeg <- Sys.which("ffmpeg")
if (nzchar(ffmpeg)) {
  cmd <- sprintf('"%s" -y -framerate %d -i "%s/frame-%%06d.png" -c:v libx264 -pix_fmt yuv420p "%s"',
                 ffmpeg, fps, out_dir, mp4_path)
  status <- system(cmd)
  if (status == 0) {
    message("MP4 written to: ", normalizePath(mp4_path, winslash = "/", mustWork = FALSE))
  } else {
    message("ffmpeg failed; PNG sequence is still available.")
  }
} else {
  message("ffmpeg not found; only PNG sequence was generated.")
}

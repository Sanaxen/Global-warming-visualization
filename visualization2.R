org_libpath <- .libPaths()

curdir = getwd()


install_libpath = paste(curdir, "/library", sep="")

.libPaths( c(install_libpath))

if (F)
{
install.packages("ncdf4", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("readr", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("R.utils", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("viridisLite", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("ragg", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("maps", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("sf", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("rnaturalearth", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
install.packages("maps", repos = "http://cran.us.r-project.org",dependencies=TRUE, lib=install_libpath, type = "binary")
}

.libPaths( org_libpath)
library(ncdf4)
library(ggplot2)
library(R.utils)
library(viridisLite)
library(maps)

url_nc_gz <- "https://data.giss.nasa.gov/pub/gistemp/gistemp1200_GHCNv4_ERSSTv5.nc.gz"

dest_gz <- file.path(tempdir(), "gistemp1200_GHCNv4_ERSSTv5.nc.gz")
dest_nc <- sub("\\.gz$", "", dest_gz)

if (!file.exists(dest_nc)) {
  download.file(url_nc_gz, destfile = dest_gz, mode = "wb")
  R.utils::gunzip(dest_gz, destname = dest_nc, overwrite = TRUE)
}

nc <- ncdf4::nc_open(dest_nc)

var_names <- names(nc$var)
dim_names <- names(nc$dim)

cat("Variables:\n")
print(var_names)
cat("Dimensions:\n")
print(dim_names)

lon_name  <- intersect(c("lon", "longitude", "X"), dim_names)
lat_name  <- intersect(c("lat", "latitude", "Y"), dim_names)
time_name <- intersect(c("time", "month", "T"), dim_names)

if (length(lon_name) == 0)  stop("lon not found")
if (length(lat_name) == 0)  stop("lat not found")
if (length(time_name) == 0) stop("time not found")

lon_name  <- lon_name[1]
lat_name  <- lat_name[1]
time_name <- time_name[1]

lon  <- ncvar_get(nc, lon_name)
lat  <- ncvar_get(nc, lat_name)
time <- ncvar_get(nc, time_name)

candidate_vars <- var_names[sapply(var_names, function(v) {
  vd <- nc$var[[v]]
  length(vd$dim) >= 3
})]

if (length(candidate_vars) == 0) {
  stop("3dim not found")
}

cat("Candidate data variables:\n")
print(candidate_vars)

varname <- candidate_vars[1]
cat("Using variable:", varname, "\n")

arr <- ncvar_get(nc, varname)
arr_dim <- dim(arr)
cat("Array dim:\n")
print(arr_dim)

vdim_names <- sapply(nc$var[[varname]]$dim, function(d) d$name)
print(vdim_names)

idx_lon  <- match(lon_name,  vdim_names)
idx_lat  <- match(lat_name,  vdim_names)
idx_time <- match(time_name, vdim_names)

if (any(is.na(c(idx_lon, idx_lat, idx_time)))) {
  stop("lon/lat/time not found")
}

perm <- c(idx_lon, idx_lat, idx_time, setdiff(seq_along(arr_dim), c(idx_lon, idx_lat, idx_time)))
arr2 <- aperm(arr, perm)

while (length(dim(arr2)) > 3) {
  arr2 <- arr2[, , , 1, drop = TRUE]
}

if (!identical(dim(arr2)[1], length(lon))) {
  stop("lon dimensions do not match.")
}
if (!identical(dim(arr2)[2], length(lat))) {
  stop("lat dimensions do not match.")
}
if (!identical(dim(arr2)[3], length(time))) {
  stop("time dimensions do not match.")
}

fill_value <- ncatt_get(nc, varname, "_FillValue")$value
missing_value <- ncatt_get(nc, varname, "missing_value")$value

if (!is.null(fill_value) && !all(is.na(fill_value))) {
  arr2[arr2 %in% fill_value] <- NA
}
if (!is.null(missing_value) && !all(is.na(missing_value))) {
  arr2[arr2 %in% missing_value] <- NA
}

parse_nc_time <- function(time_vals, time_units) {
  if (is.null(time_units) || is.na(time_units) || !nzchar(time_units)) {
    return(rep(NA_character_, length(time_vals)))
  }

  tu <- tolower(time_units)

  if (grepl("since", tu)) {
    origin_txt <- sub(".*since\\s+", "", tu)
    origin_txt <- trimws(origin_txt)

    if (grepl("^day", tu)) {
      origin_date <- as.POSIXct(origin_txt, tz = "UTC")
      if (is.na(origin_date)) origin_date <- as.POSIXct(substr(origin_txt, 1, 10), tz = "UTC")
      tt <- origin_date + time_vals * 86400
      return(format(tt, "%Y-%m"))
    }

    if (grepl("^hour", tu)) {
      origin_date <- as.POSIXct(origin_txt, tz = "UTC")
      if (is.na(origin_date)) origin_date <- as.POSIXct(substr(origin_txt, 1, 10), tz = "UTC")
      tt <- origin_date + time_vals * 3600
      return(format(tt, "%Y-%m"))
    }

    if (grepl("^month", tu)) {
      origin_date <- as.Date(substr(origin_txt, 1, 10))
      y0 <- as.integer(format(origin_date, "%Y"))
      m0 <- as.integer(format(origin_date, "%m"))

      out <- character(length(time_vals))
      for (i in seq_along(time_vals)) {
        n <- round(time_vals[i])
        y <- y0 + (m0 - 1 + n) %/% 12
        m <- (m0 - 1 + n) %% 12 + 1
        out[i] <- sprintf("%04d-%02d", y, m)
      }
      return(out)
    }
  }

  rep(NA_character_, length(time_vals))
}

time_units <- ncatt_get(nc, time_name, "units")$value
time_labels <- parse_nc_time(time, time_units)

if (all(is.na(time_labels))) {
  time_labels <- sprintf("index-%06d", seq_along(time))
}

world_map <- map_data("world")

# Move lon to [-180, 180] in case the line around the International Date Line looks strange.
normalize_lon <- function(x) {
  ((x + 180) %% 360) - 180
}

lon2 <- normalize_lon(lon)
ord_lon <- order(lon2)
lon2 <- lon2[ord_lon]
arr2 <- arr2[ord_lon, , , drop = FALSE]

if (length(lat) >= 2 && lat[1] > lat[length(lat)]) {
  lat <- rev(lat)
  arr2 <- arr2[, ncol(arr2[, , 1]):1, , drop = FALSE]
}

make_df_for_k <- function(k, arr3, lon, lat) {
  mat <- arr3[, , k]

  df <- expand.grid(
    lon = lon,
    lat = lat
  )
  df$anomaly <- as.vector(mat)
  df
}

plot_one_frame <- function(df, title_txt = NULL, zlim = c(-4, 4)) {
  ggplot(df, aes(x = lon, y = lat, fill = anomaly)) +
    geom_raster() +
    geom_path(
      data = world_map,
      aes(x = long, y = lat, group = group),
      inherit.aes = FALSE,
      linewidth = 0.22,
      color = "white",
      lineend = "round"
    ) +
    coord_quickmap(
      xlim = c(-180, 180),
      ylim = range(df$lat, na.rm = TRUE),
      expand = FALSE
    ) +
    scale_fill_gradientn(
      colours = c(
        "#08306b", "#2171b5", "#6baed6", "#c6dbef",
        "#f7f7f7",
        "#fcbba1", "#fb6a4a", "#de2d26", "#67000d"
      ),
      limits = zlim,
      oob = scales::squish,
      na.value = "grey80",
      name = "°C"
    ) +
    labs(
      title = title_txt,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

out_dir <- file.path(getwd(), "gistemp_frames_with_coast")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

zlim <- c(-4, 4)
use_ragg <- requireNamespace("ragg", quietly = TRUE)

n_time <- dim(arr2)[3]
cat("Saving", n_time, "frames to:", out_dir, "\n")

for (k in seq_len(n_time)) {
  df <- make_df_for_k(k, arr2, lon2, lat)

  p <- plot_one_frame(
    df,
    title_txt = paste("NASA GISTEMP anomaly", time_labels[k]),
    zlim = zlim
  )

  outfile <- file.path(out_dir, sprintf("frame-%06d.png", k))

  if (use_ragg) {
    ragg::agg_png(
      filename = outfile,
      width = 1600,
      height = 900,
      units = "px",
      res = 144
    )
    print(p)
    dev.off()
  } else {
    png(
      filename = outfile,
      width = 1600,
      height = 900,
      res = 144
    )
    print(p)
    dev.off()
  }

  if (k %% 50 == 0 || k == n_time) {
    cat("saved:", k, "/", n_time, "\n")
  }
}

ncdf4::nc_close(nc)

cat("Done.\n")
cat("Output folder:\n")
cat(out_dir, "\n")
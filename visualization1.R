org_libpath <- .libPaths()

curdir = getwd()


install_libpath = paste(curdir, "/library", sep="")

.libPaths( c(install_libpath))

if (T)
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

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(ncdf4)
library(R.utils)
library(viridisLite)
library(maps)
library(sf)
library(rnaturalearth)

url_global <- "https://data.giss.nasa.gov/gistemp/tabledata_v4/GLB.Ts+dSST.csv"

raw <- read_csv(
  url_global,
  skip = 1,
  show_col_types = FALSE,
  na = c("***", "****", "*****", "******", "*******")
)

print(names(raw))
print(head(raw, 3))

# Year, Jan, Feb, ..., Dec, J-D, D-N, DJF, MAM, JJA, SON
monthly <- raw %>%
  select(Year, Jan:Dec) %>%
  pivot_longer(
    cols = Jan:Dec,
    names_to = "month_abbr",
    values_to = "anomaly"
  ) %>%
  mutate(
    month = match(month_abbr, month.abb),
    date = as.Date(sprintf("%04d-%02d-01", Year, month))
  ) %>%
  arrange(date)

annual <- raw %>%
  transmute(
    year = Year,
    annual_anomaly = `J-D`
  )

print(head(monthly, 15))
print(tail(monthly, 15))
print(tail(annual, 15))

p1<-ggplot(monthly, aes(date, anomaly)) +
  geom_line() +
  labs(
    title = "NASA GISTEMP Global Temperature Anomaly",
    subtitle = "Baseline: 1951-1980",
    x = NULL,
    y = "Temperature anomaly (°C)"
  ) +
  theme_minimal(base_size = 13)

ragg::agg_png(
  filename = "NASA-GISTEMP-Global-Temperature-Anomaly1.png",
  width = 1600,
  height = 900,
  units = "px",
  res = 144
)
print(p1)
dev.off()


p2 <- ggplot(annual, aes(year, annual_anomaly)) +
  geom_col() +
  labs(
    title = "Annual Global Temperature Anomaly",
    subtitle = "NASA GISTEMP v4, baseline 1951-1980",
    x = NULL,
    y = "Temperature anomaly (°C)"
  ) +
  theme_minimal(base_size = 13)

ragg::agg_png(
  filename = "NASA-GISTEMP-Global-Temperature-Anomaly2.png",
  width = 1600,
  height = 900,
  units = "px",
  res = 144
)
print(p2)
dev.off()



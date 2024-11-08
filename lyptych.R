library(terra)

options(parallelly.fork.enable = TRUE, future.rng.onMisuse = "ignore")
library(furrr)


topo <- project(rast(sds::gebco()), rast(), by_util = TRUE)
cells <- sample(which(values(topo)[,1] > 0), 1000)

## build a simple set of 1-degree extents from a world grid
p2ex <- function(x) {
  x <- floor(x)
  cbind(x[,1], x[,1] + 1, x[,2], x[,2] + 1)
  #rep(floor(x), 2) + c(0, 0, 1, 1)
}
extents <- xyFromCell(topo, cells) |> p2ex() #|> vaster::plot_extent(add = T)

l <- vector("list")

get_stac <- function(ex) {
  h <- try(tacky:::hrefs(sds::stacit(ex, "2024-11")), silent = TRUE)
  if (!inherits(h, "try-error") && nrow(h) > 0) {
    #h$cell <- cell
    return(h)
  }
  return(NULL)
}
plan(multicore)
l <- future_map_dfr(split(t(extents), rep(1:nrow(extents), each = 4)), 
                    get_stac)



## bbox to wk rct geometry
bb2rct <- function(x) {
  exy <- do.call(rbind, x)
  wch <- exy[,3] < exy[,1]
  if (any(wch)) {
    exy[wch, 3] <-  exy[wch, 3] + 360
  }
  wk::rct(exy[,1], exy[,2], exy[,3], exy[,4])
}
x <- bb2rct(l$bbox)




library(wk)
library(geos)
rc <- as_rct(grd(rct(-180, -90, 180, 90), dx = 1, dy = 1, type = "polygons"))

tree <- geos_basic_strtree(rc)
idx <- geos_basic_strtree_query(tree, x)

m <- maps::map(plot = F)


ql <- function (x, dim = NULL, ...) 
{
  x <- x[1]
  if (!grepl("<VRT", x) && grepl("^http", x)) 
    x <- sprintf("/vsicurl/%s", x)
  info <- vapour::vapour_raster_info(x)
  dim <- tail(info$overviews, 2)
  if (is.null(dim)) 
    dim <- c(1024, 0)
  d <- gdal_raster_nara(x[1], target_dim = dim, bands = 1:info$bands, target_crs = "EPSG:4326")
  dplot <- d
  if (info$bands > 3) {
    dplot <- d[1:3]
    attributes(dplot) <- attributes(d)
  }
  try(ximage::ximage(dplot, asp = 1, ...), silent = TRUE)
  invisible(d)
}

ic <- sample(1:nrow(idx), 1)
p <- wk_coords(geos::geos_centroid(x[idx[ic, 1]]))

plot(c(x[idx[ic, 1]], rc[idx[ic, 2]]), col = c("grey", "transparent"), asp = 1/cos(p[,"y"] * pi/180))
ql(dsn::vsicurl(l$visual[idx[ic, 1]]), add = TRUE)

abline(v = 180, lty = 2)
maps::map(m, add = TRUE)


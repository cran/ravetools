## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 3.2,
  dpi = 96,
  out.width = "100%"
)
library(ravetools)
oldpar <- graphics::par(no.readonly = TRUE)

## ----primitives---------------------------------------------------------------
sphere <- vcg_sphere(sub_division = 4)
sphere

plane <- plane_geometry(width = 2.6, height = 2.6, shape = c(40, 40))
plane

## ----tilt---------------------------------------------------------------------
tilt <- new_matrix4()$make_rotation_y(25 * pi / 180)$to_array()
plane$vb <- (tilt %*% rbind(plane$vb[1:3, ], 1))[1:3, ] + c(0, 0, 0.35)

## ----normals------------------------------------------------------------------
plane <- vcg_update_normals(plane, weight = "area")
dim(plane$normals)

## ----collision----------------------------------------------------------------
cut <- vcg_detect_collision(sphere, plane, radius = 0.02)
cut$summary$y$unit_type
sum(cut$hit_unit)
head(cut$representation)

## ----band---------------------------------------------------------------------
band <- plane
band$it <- plane$it[, which(cut$hit_unit), drop = FALSE]
band

## ----radius-------------------------------------------------------------------
vapply(c(0.02, 0.1, 0.3), function(r) {
  sum(vcg_detect_collision(sphere, plane, radius = r)$hit_unit)
}, 0L)

## ----streamlines--------------------------------------------------------------
streamlines <- rbind(
  cbind(seq(-3, 3, by = 0.5), 0, 0),
  c(NA, NA, NA),
  cbind(seq(-3, 3, by = 0.5), 5, 0)
)
res <- vcg_detect_collision(sphere, streamlines,
                            mode_y = "segments", radius = 0.1)
res$hit_unit
res$representation

## ----test-level---------------------------------------------------------------
vcg_detect_collision(sphere, streamlines, mode_y = "segments",
                     radius = 0.1, test_level = "whole")$collide

## ----interior-----------------------------------------------------------------
center <- rbind(c(0, 0, 0))
vcg_detect_collision(sphere, center, radius = 0.1)$collide
vcg_detect_collision(sphere, center, radius = 0.1,
                     include_interior = TRUE)$collide

## ----kdtree-------------------------------------------------------------------
kd <- vcg_kdtree_nearest(target = sphere, query = plane, k = 1)
str(kd)
range(kd$distance)

## ----kdtree-k-----------------------------------------------------------------
kd3 <- vcg_kdtree_nearest(target = sphere, query = plane, k = 3)
head(kd3$index)

## ----raycast------------------------------------------------------------------
above <- plane$vb
above[3, ] <- above[3, ] + 1.8
rays <- vcg_raycaster(sphere, ray_origin = above, ray_direction = c(0, 0, -1))
str(rays[c("has_intersection", "distance", "face_index")])
sum(rays$has_intersection)

## ----pierce-------------------------------------------------------------------
pierce <- structure(
  list(vb = rays$intersection[, rays$has_intersection, drop = FALSE]),
  class = "mesh3d"
)

## ----query-figure, fig.height = 3.6-------------------------------------------
col <- color_ramp_continuous(kd$distance[, 1],
                             cmap = c("#f2f2f2", "#7fa8c9", "#2c5f8a"))

graphics::par(mfrow = c(1, 2), mar = c(0.1, 0.1, 2.1, 0.1), cex.main = 0.95)
plot_mesh_polygon(
  list(sphere, plane, band),
  col = list("gray55", col, "#a33a3a"), alpha = c(0.5, 0.92, 1),
  eye = c(3.5, -4, 2.2), up = c(0, 0, 1), zoom = 1.15,
  shadow_color = "white", ambient_intensity = 0.55,
  main = "Distance; cut faces in red"
)
plot_mesh_polygon(
  list(sphere, pierce),
  col = list("gray72", "#a33a3a"), cex = 0.03,
  eye = c(3.5, -4, 2.2), up = c(0, 0, 1), zoom = 0.68,
  shadow_color = "white", ambient_intensity = 0.55,
  main = "Ray entry points"
)

## ----isosurface---------------------------------------------------------------
data("left_hippocampus_mask", package = "ravetools")
dim(left_hippocampus_mask)

raw_mesh <- vcg_isosurface(left_hippocampus_mask)
raw_mesh

## ----mesh-from-volume---------------------------------------------------------
smoothed <- mesh_from_volume(
  left_hippocampus_mask, output_format = "rgl", threshold = 0.5,
  remesh = TRUE, remesh_voxel_size = 1, smooth = TRUE, verbose = FALSE
)
smoothed

## ----defects------------------------------------------------------------------
vcg_count_edge_defects(raw_mesh)

## ----measure------------------------------------------------------------------
vcg_mesh_volume(raw_mesh)
vcg_average_edge_length(raw_mesh)
vcg_max_edge_length(raw_mesh)

## ----curvature-fields---------------------------------------------------------
curv <- mris_curvature(raw_mesh)
str(curv)

## ----fix----------------------------------------------------------------------
mesh <- vcg_fix_defects(raw_mesh, verbose = FALSE)
info <- attr(mesh, "info")
info[c("boundary_edges_before", "boundary_edges_after",
       "holes_filled", "is_closed_manifold")]

# the same call that warned above, now on a closed surface
vcg_mesh_volume(mesh)

## ----center-------------------------------------------------------------------
mesh$vb[1:3, ] <- mesh$vb[1:3, ] - rowMeans(mesh$vb[1:3, ])

## ----remesh-------------------------------------------------------------------
uniform  <- vcg_uniform_remesh(mesh, voxel_size = 1, verbose = FALSE)
split    <- vcg_subdivision(mesh, method = "edge")
capped   <- vcg_subdivide_max_edge_length(mesh, max_edge_len = 0.8)
isotropic <- mris_remesh(mesh, target_edge_length = 1.5, verbose = FALSE)

data.frame(
  method = c("input", "vcg_uniform_remesh", "vcg_subdivision",
             "vcg_subdivide_max_edge_length", "mris_remesh"),
  vertices = c(ncol(mesh$vb), ncol(uniform$vb), ncol(split$vb),
               ncol(capped$vb), ncol(isotropic$vb)),
  avg_edge = round(vapply(list(mesh, uniform, split, capped, isotropic),
                          vcg_average_edge_length, 0), 3),
  max_edge = round(vapply(list(mesh, uniform, split, capped, isotropic),
                          vcg_max_edge_length, 0), 3)
)

## ----smooth-------------------------------------------------------------------
taubin <- vcg_smooth_explicit(mesh, type = "taubin", iteration = 10)
implicit <- vcg_smooth_implicit(mesh, lambda = 0.2, degree = 2)
fs_style <- mris_smooth(mesh, niterations = 20L)

vapply(list(mesh, taubin, implicit, fs_style), vcg_mesh_volume, 0)

## ----inflate------------------------------------------------------------------
inflated <- mris_inflate(fs_style, n_averages = 4L, niterations = 8L,
                         scale_brain = FALSE, verbose = FALSE)
names(inflated)
range(inflated$sulc)

## ----sphere-map---------------------------------------------------------------
spherical <- mris_sphere(fs_style, target_radius = 100, verbose = FALSE)
radius <- sqrt(colSums(spherical$vb[1:3, ]^2))
c(min = min(radius), max = max(radius), cv = stats::sd(radius) / mean(radius))

## ----inflate-figure, fig.height = 3.2-----------------------------------------
curv <- mris_curvature(fs_style)
lim <- stats::quantile(abs(curv$mean), 0.95)
curv_col <- color_ramp_continuous(
  curv$mean, clim = c(-lim, lim),
  cmap = c("#2c5f8a", "#f2f2f2", "#a33a3a")
)

graphics::par(mfrow = c(1, 2), mar = c(0.1, 0.1, 2.1, 0.1))
plot_mesh_polygon(fs_style, col = curv_col, eye = c(0, 100, 30),
                  up = c(0, 0, 1), zoom = 1.1, main = "Mean curvature")
plot_mesh_polygon(inflated$mesh, col = curv_col, eye = c(0, 100, 30),
                  up = c(0, 0, 1), zoom = 1.1, main = "The same, inflated")

## ----subset-------------------------------------------------------------------
selector <- mesh$vb[1, ] > 0
half <- vcg_subset_vertex(mesh, selector)
c(input = ncol(mesh$vb), kept = ncol(half$vb))

## ----patch--------------------------------------------------------------------
target <- vcg_uniform_remesh(vcg_sphere(), verbose = FALSE)
patches <- vcg_mesh_patch(target, waypoints = diag(1, 3))
vapply(patches, function(p) ncol(p$it), 0L)

## ----patch-figure, fig.height = 3.4-------------------------------------------
graphics::par(mar = c(0.1, 0.1, 0.1, 0.1))
plot_mesh_polygon(patches, col = list("#a33a3a", "gray70"),
                  alpha = c(1, 0.55), eye = c(10, 10, 10), zoom = 1.2,
                  shadow_color = "white", ambient_intensity = 0.55)

## ----dijkstra-----------------------------------------------------------------
dist <- dijkstras_surface_distance(
  positions = t(mesh$vb[1:3, ]),
  faces = t(mesh$it),
  start_node = 1,
  face_index_start = 1
)
path <- surface_path(dist, target_node = ncol(mesh$vb))
c(vertices_on_path = length(path$path), length = max(path$distance))

## ----clipping, fig.height = 3.4-----------------------------------------------
bbox <- apply(fs_style$vb[1:3, ], 1L, range)
set.seed(1)
candidates <- cbind(
  stats::runif(600, bbox[1, 1], bbox[2, 1]),
  stats::runif(600, bbox[1, 2], bbox[2, 2]),
  stats::runif(600, bbox[1, 3], bbox[2, 3])
)
inside <- vcg_detect_collision(fs_style, candidates,
                               include_interior = TRUE)$hit_unit %in% TRUE
probes <- structure(
  list(vb = t(candidates[inside, , drop = FALSE])),
  class = "mesh3d"
)
sum(inside)

eye <- c(0, 100, 30)
graphics::par(mfrow = c(1, 3), mar = c(0.1, 0.1, 2.1, 0.1), cex.main = 0.95)

plot_mesh_polygon(fs_style, col = "steelblue", eye = eye, up = c(0, 0, 1),
                  zoom = 1.1, main = "Whole surface")

plot_mesh_polygon(
  list(fs_style, probes),
  col = list("steelblue", "#a33a3a"), cex = 1.1,
  eye = eye, up = c(0, 0, 1), zoom = 1.1,
  clipping_plane = c(0, 1, 0, 0, -1),
  clipping_plane_enabled = c(TRUE, FALSE),
  main = "Cut open, probes exempt"
)

# `plot_mesh_dotcloud` has no `main`; add the title afterwards
plot_mesh_dotcloud(fs_style, col = "steelblue", eye = eye, up = c(0, 0, 1),
                   zoom = 1.1, cex = 0.45)
graphics::title(main = "plot_mesh_dotcloud")

## ----rgl, eval = FALSE--------------------------------------------------------
# rgl_view({
#   rgl_call("shade3d", mesh, col = "steelblue")
#   rgl_call("wire3d", mesh, col = "black")
# })

## ----cleanup, include = FALSE-------------------------------------------------
graphics::par(oldpar)


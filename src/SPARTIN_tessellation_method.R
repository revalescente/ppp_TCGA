library(SPARTIN)

library(spatstat)
set.seed(10000)

CustomHeatmap = CH = function(m, t="", min.v = NA, max.v = NA){
  plot.tib = tibble(r = numeric(0), c = numeric(0), val = numeric(0))
  for(i in 1:nrow(m)){
    plot.tib = bind_rows(plot.tib, tibble(
      r = rep(i, ncol(m)),
      c = 1:ncol(m),
      val = as.numeric(m[i,])))
  }
  
  plot.tib$c = as.factor(plot.tib$c)
  
  if(is.na(min.v) || is.na(max.v)){
    min.v = min(plot.tib$val, na.rm = T)
    max.v = max(plot.tib$val, na.rm = T)
  }
  cols = colorRampPalette(c("#0a0722", "#3d0965", "#721a6e", "#a52c60", "#d44842",
                            "#f37819", "#fcb216"))(7) #, "#f1f179"))(8)
  
  ggplot(plot.tib) +
    geom_tile(mapping = aes(x = c, y = reorder(r, -r), fill = val)) +
    labs(fill = "") +
    ggtitle(t) +
    scale_fill_gradientn(colors = cols,
                         limits = c(min.v, max.v),
                         # breaks = seq(min.v, max.v, 2),
                         na.value = 'white') +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.background = element_rect(fill = "white",
                                          colour = "black"),
          panel.border = element_rect(linetype = "solid", fill = NA))
}

# Visualize PPP
vis = function(p, t="", suppWarn = F){
  if(suppWarn){
    suppressWarnings(
      plot(p, cols = c("black", "red"),
           shape = c("circles", "circles"),
           size = 6,
           main = t)
    )
  }else{
    plot(p, cols = c("black", "red"),
         shape = c("circles", "circles"),
         size = 6,
         main = t)
  }
}

ex_ppp = SimulateData(n1 = 100, n2 = 40, phi = 0.3,
                      winX = 300, winY = 300, r = 30)
ex_ppp
plot(ex_ppp)
ex_freq = FitHSFreq(ex_ppp, r = 30, quad.spacing = 3)
coef(summary(ex_freq))

ex_bayes = FitHSBayes(ex_ppp, r = 30, quad.spacing = 3)
ex_bayes

# Estimating CTIP

ex_CTIP = CTIP(ex_ppp, r = 30, quad.spacing = 3, 
               n.null = 5, n.burn = 1000,
               n.sample = 11000, null.n.burn = 1000, null.n.sample = 11000,
               n.thin = 5)
ex_CTIP

# Tesselation and Intensity Thresholding

# Simulate a biopsy
ex_biopsy_part_1 = SimulateData(n1 = 5000, n2 = 1000, phi = 0.4,
                                winX = 1000, winY = 1000, r = 30)

ex_biopsy_part_2 = SimulateData(n1 = 5000, n2 = 1000, phi = 0.4,
                                winX = 1000, winY = 1000, r = 30)

# Place in a large window to demonstrate intensity thresholding in action
ex_biopsy_excess = ppp(x = c(ex_biopsy_part_1$x, ex_biopsy_part_2$x - 1000),
                       y = c(ex_biopsy_part_1$y, ex_biopsy_part_2$y - 1000),
                       marks = c(ex_biopsy_part_1$marks,
                                 ex_biopsy_part_2$marks),
                       window = owin(xrange = c(-1000, 1000),
                                     yrange = c(-1000, 1000)))

vis(ex_biopsy_excess)

# Bad tessellation
bad_tessellation_1 = TessellateBiopsy(PPPToTibble(ex_biopsy_excess),
                                      sigma = 1, eps = 15,
                                      threshold = (ex_biopsy_excess$n)/
                                        area.owin(ex_biopsy_excess$window),
                                      clust.size = 75,
                                      max.clust.size = 100)
vis(bad_tessellation_1$window)
vis(bad_tessellation_1$tiles[[1]])

# Good tessellation
good_tessellation = TessellateBiopsy(PPPToTibble(ex_biopsy_excess),
                                     sigma = 30, eps = 15,
                                     threshold = (ex_biopsy_excess$n)/
                                       area.owin(ex_biopsy_excess$window),
                                     clust.size = 75,
                                     max.clust.size = 100)
vis(good_tessellation$window)
vis(good_tessellation$tiles[[2]])

#visualization of all the tiles
bad_tes_n_tum = purrr::map_dbl(bad_tessellation_1$tiles, ~ sum(.x$marks == 1))
good_tes_n_tum = purrr::map_dbl(good_tessellation$tiles, ~ sum(.x$marks == 1))

PlotTessellation(bad_tessellation_1, bad_tes_n_tum, "# Tumor Cells")
PlotTessellation(good_tessellation, good_tes_n_tum, "# Tumor Cells")

# WSI data
library(ggplot2)
library(dplyr)
data("example_WSI")
head(example_WSI)


ggplot(example_WSI |> filter(Mark == 1)) + 
  geom_point(aes(x = CentroidX, y = CentroidY), size = 0.1, color = "black") + 
  ggtitle("Tumor cell locations") + 
  theme_bw()

ggplot(example_WSI |> filter(Mark == 2)) + 
  geom_point(aes(x = CentroidX, y = CentroidY), size = 0.1, color = "red") + 
  ggtitle("Immune cell locations") + 
  theme_bw()

left_chunk = example_WSI %>% 
  filter(CentroidX < 10000)

ggplot(left_chunk %>% filter(Mark == 1)) + 
  geom_point(aes(x = CentroidX, y = CentroidY, color = Mark), size = 0.1, color = "black") + 
  ggtitle("Tumor cell locations") + 
  theme_bw()

ggplot(left_chunk %>% filter(Mark == 2)) + 
  geom_point(aes(x = CentroidX, y = CentroidY, color = Mark), size = 0.1, color = "red") + 
  ggtitle("Immune cell locations") + 
  theme_bw()

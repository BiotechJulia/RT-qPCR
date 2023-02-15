# Set the working directory (this is where your source files are and where the output figures will be saved).
setwd(choose.dir(getwd()))

# Install these packages, if you don't have them installed.
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("reshape2")
install.packages("dplyr")
install.packages("remotes")
remotes::install_github("wilkelab/ggtext")
install.packages("ggpubr")
install.packages("multcompView")
install.packages("ggthemes")

# Load libraries.
library(ggplot2)
library(tidyverse)
library(reshape2)
library(dplyr)
library(remotes)
library(ggtext)
library(ggpubr)
library(multcompView)
library(ggthemes)

# Load RT-qPCR data prepared accordingly to the sample data. 
# The values here are mean CT values from 3 technical replicates. 3 biological replicates were tested.
dane <- read.csv("sample data.csv",sep = ",")

# Calculate target gene expression relatively to reference gene expression. 
# Here, the name of reference gene is "REF" and the name of target gene is "TARGET1". These names are used in the code.
# If you want to set other names, change them in the code too.
ref <- subset(dane, dane$Target.Name=="REF")
target1 <- subset(dane, dane$Target.Name=="TARGET1")
target2 <- subset(dane, dane$Target.Name=="TARGET2")
rel_expr <- 2^ref$C./2^target1$C.
dane1 <- cbind(target1, rel_expr)
rel_expr <- 2^ref$C./2^target2$C.
dane2 <- cbind(target2, rel_expr)

# Now do the statistics - here I will show you how to perform ANOVA with TukeyHSD post-hoc.
# TARGET1
anova <- aov(rel_expr ~ Sample.Name*Condition, data = dane1)
tukey <- TukeyHSD(anova)
cld1 <- multcompLetters4(anova, tukey)
cld1 <- as.data.frame.list(cld1$`Sample.Name:Condition`)
tukey$`Sample.Name:Condition`

# TARGET2
anova <- aov(rel_expr ~ Sample.Name*Condition, data = dane2)
tukey <- TukeyHSD(anova)
cld2 <- multcompLetters4(anova, tukey)
cld2 <- as.data.frame.list(cld2$`Sample.Name:Condition`)
tukey$`Sample.Name:Condition`

# Prepare the data to plot - calculate mean and standard deviation for obtained relative expression values.
# TARGET1
new1 <- dane1 %>% group_by(Condition, Sample.Name, Target.Name) %>% summarise_each(funs(mean, sd))

new1$lower <- with(new1, rel_expr_mean-rel_expr_sd)
new1$upper <- with(new1, rel_expr_mean+rel_expr_sd)
new1$sample <- paste(new1$Sample.Name,":",new1$Condition, sep="")
cld1 <- cld1[ order((row.names(cld1))), ,drop=F]
new1 <- new1[ order((new1$sample)), ,drop=F]
new1$cld <- cld1$Letters

# TARGET2
new2 <- dane2 %>% group_by(Condition, Sample.Name, Target.Name) %>% summarise_each(funs(mean, sd))

new2$lower <- with(new2, rel_expr_mean-rel_expr_sd)
new2$upper <- with(new2, rel_expr_mean+rel_expr_sd)
new2$sample <- paste(new2$Sample.Name,":",new2$Condition, sep="")
cld2 <- cld2[ order((row.names(cld2))), ,drop=F]
new2 <- new2[ order((new2$sample)), ,drop=F]
new2$cld <- cld2$Letters

# If you want, change the order of samples that will show on X axis of the plot and/or assign new labels for them.
# In the same way you can reorder conditions.
new <- rbind(new1, new2)
new$Sample.Name <- factor (new$Sample.Name, levels = c("wt","mutant"), labels = expression("WT","*mutant*"))

######################################## PLOT ###################################################
stack <- ggplot(new, aes(x=Sample.Name, y=rel_expr_mean, fill=Condition)) +
  facet_wrap(~new$Target.Name, scales = "free") +
  scale_fill_manual(values = c("white",  "black")) +
  geom_bar(stat = "identity",colour="black",width=.7,size=.5, position = position_dodge(.7)) +
  scale_y_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 1), limits = c(0,(max(new$rel_expr_mean)*1.5))) +
  theme_classic() +
  labs(
    title=" ",
    x=" ",
    y="relative expression"
  )  + 
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.1,size=.5, position = position_dodge(.7)) + 
  theme(legend.title = element_blank(), 
        axis.text.y = element_text(color="black"), 
        axis.text.x = element_text(color="black")) +
  geom_text(aes(label = cld, y = rel_expr_mean+rel_expr_sd),size = 6*0.352777778, vjust = -0.5, position=position_dodge(.7)) +
  # scale_x_discrete(expand = c(0.2, 0)) + 
  theme(text = element_text(size = 6),strip.text.x = element_text(
    size = 6, color = "black", face = "italic"
  ),
  axis.text.x = element_markdown(), axis.ticks = (element_line(colour="black"))) +
  theme(legend.key.size = unit(0.1, 'in'), plot.title = element_text(size = 8), legend.position = c(0.12, 0.86),  legend.background=element_rect(fill = alpha("white", 0)))+ 
  theme(axis.text.x = element_markdown(angle = 0), plot.margin = margin(0, 0, 0, 0, "cm"))


stack

# As you can see, two target genes differ in the expression levels.
# This is a function that lets you set separate scales for each plot
scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)


facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) ||
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}

# Lets draw the plots again with separate scales
a <- stack + facet_wrap_custom(~Target.Name, scales = "free", ncol = 2, scale_overrides = list(
  scale_override(1, scale_y_continuous(limits = c(0, 300), expand = c(0,0))),
  scale_override(2, scale_y_continuous(limits = c(0, 400), expand = c(0,0)))
))

a

# Save the plot.
# .svg format file lets you edit it with vector graphic editors like Inkscape
svg(filename = "Test RT-qPCR.svg", width = 2, height = 1.5)
a
dev.off()

# .jpg
jpeg("Test RT-qPCR.jpeg", width = 2,
     height    = 1.5,
     units     = "in",
     res       = 600,
     pointsize = 6)
a
dev.off()

# End of the code.
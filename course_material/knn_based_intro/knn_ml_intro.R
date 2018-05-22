## ----style, echo=FALSE, results="asis", cache=FALSE----------------------
library("knitr")
options(digits = 2, width = 80)
golden_ratio <- (1 + sqrt(5)) / 2
opts_chunk$set(echo = TRUE, tidy = FALSE, include = TRUE,
               dev=c('png', 'pdf', 'svg'), fig.height = 5, fig.width = 4 * golden_ratio, comment = '  ', dpi = 300,
cache = TRUE)

## ---- echo=FALSE, cache=FALSE--------------------------------------------
print(date())

## ----setup, cache = FALSE------------------------------------------------
library("rmarkdown")
library("BiocStyle")
library("magrittr")
library("stringr")
library("ggthemes")
library("scales")
library("ggbeeswarm")
library("factoextra")
library("tidyverse")
library("readxl")
library("ggrepel")
library("RColorBrewer")
library("mlr")
library("psych")

set.seed(222)

theme_set(theme_solarized(base_size = 18))

## ----importAnnotation----------------------------------------------------
plate_map <- read_excel("plate_mapping.xlsx")
plate_map

## ----importDataTable, eval=TRUE------------------------------------------
load("raw_data.RData")

## ----reshape, dependson="import_data_table"------------------------------
tidy_raw_data  <- rownames_to_column(as.data.frame(raw_data), 
                                     var = "class") %>%
                  gather(key = "well", value = "count", -class)

tidy_raw_data

tidy_raw_data$well <- str_replace(tidy_raw_data$well, "^W([A-H][0-9]{2})_P1", "\\1_01")

#join annotation

input_data <- left_join(tidy_raw_data, plate_map, by = c("well" = "Position"))

input_data

## ----groupingAndSummarizing, dependson="reshape"-------------------------

no_cells_per_well <- input_data %>%
                     group_by(well) %>%
                     dplyr::summarize(no_cells = sum(count))

no_cells_per_well

data_with_sums <-  left_join(input_data, no_cells_per_well)

data_processed <- mutate(data_with_sums, perc = count / no_cells, 
                       z_score = logit(perc))

data_processed

## ----iputDataAfterPreP---------------------------------------------------
data_processed

## ----getDataForML--------------------------------------------------------
data_for_ML <- select(data_processed, class, Well, Group, z_score) %>%
               spread(class, z_score)
rownames(data_for_ML) <- data_for_ML$Well
data_for_ML$Group %<>% factor(levels = c("neg", "target"))
data_for_ML$Well <- NULL

## ----addClassColumn------------------------------------------------------
data_for_ML

ggplot(data = data_for_ML, aes(x = Group, y = inter, color = Group)) +
      geom_beeswarm(cex = 3) +
      ggtitle("Interphase z-scores") +
      scale_color_brewer(palette = "Accent")

## ----PCA, dependson="grouping_and_summarizing", fig.cap="PCA plot"-------
data_for_PCA <- data_processed %>% 
                dplyr::select(class, well, z_score) %>%
                spread(key = class, value = z_score)

PCA <- prcomp(data_for_PCA[, -1], center = TRUE, scale. = TRUE)


genes <- input_data %>%
         group_by(well) %>%
          dplyr::summarize(gene = unique(`Gene Symbol`))

genes <- ifelse(is.na(genes$gene), "empty", genes$gene)

dataGG = data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                    PC3 = PCA$x[,3], PC4 = PCA$x[,4],
                    genes)

pl <- (ggplot(dataGG, aes(x = PC1, y = PC2, color =  genes))
      + geom_text(aes(label = genes), size = I(2))
      + coord_fixed(ratio = (PCA$sdev)[2] / (PCA$sdev)[1])
      + ggtitle("Principal components plot")
      )

pl + scale_color_tableau(palette = "tableau20")

## ----plotly, eval=FALSE, fig.width= 4.5, dependson="PCA"-----------------
   if(!("plotly" %in% installed.packages()[,1])){
     install.packages("plotly")
   }
   
   library(plotly)
   ggplotly(pl)

## ----varImp, dependson="PCA"---------------------------------------------
loadings <- PCA$rotation[, 1:2]
loadings_gg <- loadings %>%
               as.data.frame() %>%
               rownames_to_column(var = "class") %>%
               dplyr::select(class, PC1, PC2) %>%
               gather( key = "comp", value = "loading", PC1:PC2)
  
ggplot(loadings_gg, aes(x = class, y = loading, fill = class)) +
      facet_wrap( ~ comp) +
      geom_bar(stat = "identity", position = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_fill_brewer(type = "qual", palette = "Set3") 


## ----coCircle, dependson="var_imp", fig.small=TRUE, fig.cap="Correlation Circle"----
fviz_pca_var(PCA, col.circle = "black", 
             title = "Correlation Circle for the PCA") + 
              coord_equal() 

## ----findKNN, warning=FALSE----------------------------------------------
cl_task <- makeClassifTask(id = "HTScreen", data = data_for_ML,
                           target = "Group", positive = "target")
cl_task

all_cl <- listLearners("classif")
filter(all_cl, str_detect("knn", all_cl$short.name))

## ----setHyper------------------------------------------------------------
lrn <- makeLearner("classif.knn", k = 5)
print(lrn)
getHyperPars(lrn) 

## ----trainEvaluate-------------------------------------------------------

n <- getTaskSize(cl_task)

train_set <- sample(n, size = round(n * 0.8))
head(train_set)

screen_knn <- train(lrn, cl_task, subset = train_set)
screen_knn

test_set <- setdiff(seq_len(n), train_set)
pred <- predict(screen_knn, task = cl_task, subset = test_set)

pred

performance(pred)

## ----calcCF--------------------------------------------------------------
calculateConfusionMatrix(pred = pred)

## ----vizPred, fig.cap="Vizualization of the KNN decision Boundary", fig.width= 7, fig.wide = TRUE----
learn_pl <- plotLearnerPrediction(learner = lrn, task = cl_task,
                      features = c("inter", "ana"), cv = 5) + 
  scale_fill_tableau(palette = "tableau10light") + 
  geom_tile(alpha = 0.01) +
  coord_fixed()

learn_pl 

## ---- cache=FALSE--------------------------------------------------------
sessionInfo()


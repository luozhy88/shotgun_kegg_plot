
library(tidyverse)
library(DescTools)
library(patchwork)


# width=9.0
# height=5.5
# p_value=NULL
# p_adj=0.05
# out.name="output/"
# mydata=mydata3
# zero_percent=70


# 创建sample ID
sample <- paste0("Sample", 1:100)
# 创建group
group <- sample(c("A", "B"), 100, replace = TRUE)
# 创建特征列
features <- matrix(0, nrow = 100, ncol = 8)
# 给分组A和B的特征赋予不同的均值
group_A_mean <- c(5, 6, 7, 8, 9, 10, 11, 12)
group_B_mean <- c(1, 2, 3, 4, 5, 6, 7, 8)

for (i in 1:8) {
  if (group_A_mean[i] > group_B_mean[i]) {
    features[group == "A", i] <- rnorm(sum(group == "A"), mean = group_A_mean[i])
    features[group == "B", i] <- rnorm(sum(group == "B"), mean = group_B_mean[i])
  } else {
    features[group == "A", i] <- rnorm(sum(group == "A"), mean = group_B_mean[i])
    features[group == "B", i] <- rnorm(sum(group == "B"), mean = group_A_mean[i])
  }
}

# 创建Dataframe
mydata3 <- data.frame(sample, group, features)




STAMP<- function(mydata=mydata3,zero_percent=70,p_adj=0.05,p_value=NULL,width=9.0,height=5.5,out.name="output/"){
      
    ######################################### Readme #########################################
    # 1 and 2 column are sample and group in mydata3,the order of features are ranked
    ######################################### Readme #########################################
    dir.create(out.name,recursive = T)
    ### ready for calculation
    ### remove constant
    mydata4 <- janitor::remove_constant(mydata)
    # Calculate the percentage of zero values in each column and filter columns 样本小于等于70% 含有零的特征保留
    zero_percent=70
    features_sel <- mydata4[, -c(1:2)] %>%
      summarise(across(everything(), ~ mean(. == 0) * 100)) %>%
      select(where(~ . <= zero_percent))
    
    mydata5 <- mydata4[, c("sample", "group", colnames(features_sel))]
    
    ########## ready to calculate
    colnames(mydata5) <- make.names(colnames(mydata5))
    
    mydata5 <- as.data.frame(mydata5)
    
    DF <- mydata5
    
    ### create all combinations
    pair <- DescTools::CombPairs(as.character(unique(DF$group)))
    
    for (j in 1:nrow(pair)) {
      # j=1
      DF2 <- subset(DF, group == pair[j, 1] | group == pair[j, 2])
      
      comparison <- paste0(pair[j, 1], "_", pair[j, 2])
      
      ## real mean difference
      diff.real <- DF2 %>%
        dplyr::select_if(is.numeric) %>%
        purrr::map_df(~ DescTools::MeanDiffCI(. ~ group, data = DF), .id = "var")
      
      ## model based mean difference
      diff.model <- DF2 %>%
        dplyr::select_if(is.numeric) %>%
        purrr::map_df(~ broom::tidy(wilcox.test(. ~ group, data = DF2, conf.int = T, conf.level = .95)), .id = "var")
      diff.model$p.adj <- p.adjust(diff.model$p.value, "BH")
      
      ### merge
      diff <- left_join(diff.model, diff.real)
      #    
      if(is.null(p_value) ){
        diff.sig <- dplyr::filter(diff, p.adj <= p_adj)
        # diff.sig <- subset(diff, p.adj < p.adj)
        sign_stat=glue::glue("p.adj",p.adj)
      }else{
        diff.sig <- dplyr::filter(diff, p.value <= p_value)
        sign_stat=glue::glue("p.value",p.value)
      }
      
      diff.sig$col <- ifelse(diff.sig$meandiff > 0, "01plus", "02minus")
      write.csv(diff.sig,glue::glue( out.name, sign_stat,"_Func_", comparison, ".csv"))
      ##
      DF3 <- reshape2::melt(DF2)
      DF3$variable <- as.character(DF3$variable)
      ##
      DF4 <- DF3[DF3$variable %in% diff.sig$var, ]
      
      ### take the labels by mean diff that matches the plot p2 below
      level_order <- diff.sig%>%arrange(meandiff)%>%pull(var)
      
      ## base plot
      p1 <- ggplot(data = DF4, aes(x = variable, y = value, fill = group))
      
      ### add alternated background
      for (i in 1:(nrow(diff.sig))) {
        p1 <- p1 + annotate("rect",
                            xmin = i - 0.5, xmax = i + 0.5, ymin = 0, ymax = Inf,
                            fill = ifelse(i %% 2 == 0, "white", "gray95")
        ) + scale_x_discrete()
      }
      
      p1 <- p1 +
        geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean") +
        xlab("") +
        ylab("Proportion (%)") +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black"),
          axis.ticks = element_line(color = "black", linewidth = 0.5)
        ) +
        scale_x_discrete(limits = level_order) +
        coord_flip() +
        theme(legend.position = "top")
      # Use custom colors
      p1 <- p1 + scale_fill_manual(values = c("darkorange", "deepskyblue"))
      
      
      ### base plot
      p2 <- ggplot(diff.sig, aes(x = reorder(var,+meandiff), y = meandiff, fill = col))
      ### add alternated background
      for (i in 1:(nrow(diff.sig))) {
        p2 <- p2 + annotate("rect",
                            xmin = i - 0.5, xmax = i + 0.5, ymin = -Inf, ymax = Inf,
                            fill = ifelse(i %% 2 == 0, "white", "gray95")
        ) + scale_x_discrete()
      }
      ## add the other components
      p2 <- p2 +
        geom_point(stat = "identity", size = 3, aes(col = col), position = position_dodge(width = 1)) +
        geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),
                      width = .2, linewidth = 0.5,
                      position = "dodge"
        ) +
        xlab("") +
        ylab("Difference between proportions (%)") +
        geom_hline(yintercept = 0, lty = 2) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(color = "black"),
          axis.ticks.x = element_line(color = "black", linewidth = 0.5)
        ) +
        coord_flip() +
        theme(legend.position = "none") +
        theme(axis.text.y = element_blank())
      p2 <- p2 + scale_color_manual(values = c("darkorange", "deepskyblue"))
      p1 + p2
      ##
      

      ggsave(filename = glue::glue( out.name, sign_stat,"_Func_", comparison, "_smallSize.pdf"), width = 9.0, height = height)
      ggsave(filename = glue::glue( out.name, sign_stat,"_Func_", comparison, "_mediumSize.pdf"), width = 9.0, height = height +2)
      ggsave(filename = glue::glue( out.name, sign_stat,"_Func_", comparison, "_bigSize.pdf"), width = 9.0, height = height +3)
    }
}

STAMP(mydata=mydata3,zero_percent=70,p_adj=0.05,p_value=NULL,width=9.0,height=5.5,out.name="output/")

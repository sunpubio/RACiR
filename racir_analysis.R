

# Library setting ---------------------------------------------------------
library(racir)
library(tidyverse)
library(ggplot2)
library(minpack.lm)
library(nleqslv)
library(Li6800fixer)
print("Hello")



# Folder setting ----------------------------------------------------------

getwd()

broad_dir <- "data/broad"
conifer_dir <- "data/conifer"

list_broadleaf <- list.dirs(broad_dir, full.names = T, recursive = F)
list_conifer <- list.dirs(conifer_dir, full.names = T, recursive = F)


# Conifer -----------------------------------------------------------------

for (i in list_conifer) {
  i="data/conifer/Karamatsu"
  
  all_files <- list.files(i, full.names = T, recursive = F) %>% .[!grepl("\\.xlsx?$", .)]
  aci_files <- all_files[grepl("\\_ACi", all_files)]
  area_file <- 
  
  for (m in aci_files) {
    
    m="data/conifer/Karamatsu/Karamatsu_1_ACi"
    
    ACi_data <- read_6800(m)
    plant_name <- str_extract(m, "[^/]+(?=_ACi)")
    ACi_data$PlantID <- plant_name
    
    
    # ACiデータ再計算
    area <- read.csv("leaf_area_manual_petiole_area.csv", header = TRUE) %>%
      mutate(
        # ab_数字_数字_数字 または ab_数字_ACi を抽出
        Leafname = str_extract(File, plant_name),
        str_extract()
        # Area を m² に変換
        Area_m2 = Area_mm2 / 1000000
      )
    
    
    
    
    
    
    
    names(ACi_data) <- make.unique(names(ACi_data))
    
    Judge_df <- data.frame(NULL)#Vcmax,Jmax用
    
    # 葉温決定
    t <- mean(ACi_data$Tleaf)
    
    for (j in 2:(nrow(ACi_data)-1)) {　　#点2つじゃ回帰出来なかったため、3つから行う。
      subset_data1 <- ACi_data[1:j, ] # 最初からi番目までのデータを取得
      subset_data2 <- ACi_data[(j + 1):nrow(ACi_data), ] # i+1番目から最後までのデータを取得
      print(j)
      
      
      model1 <- nlsLM(A ~ Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd, data = subset_data1, start=c(Vcmax=50, Rd=0.8), control = nls.control(maxiter = 100))
      Vcmax <- coef(model1)[1]
      Rd <- coef(model1)[2]
      
      model2 <- nlsLM(A ~ Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd, data = subset_data2, start=c(Jmax=80))
      Jmax <- coef(model2)[1]
      
      #連立方程式を解く
      Ci<- 300
      A<- 10
      initial_guess <- c(Ci, A)
      
      equations <- function(vars) {
        Ci <- vars[1]
        A <- vars[2]
        eq1 <- Vcmax * (Ci - g(t)) / (Ci + Kc(t) * (1 + O2 / Ko(t))) - Rd - A
        eq2 <- Jmax * (Ci - g(t)) / (4 * Ci + 8 * g(t)) - Rd - A
        return(c(eq1, eq2))
      }
      result <- nleqslv(initial_guess, equations)
      
      #妥当性評価
      if(all(result$x[1] <= subset_data2$Ci) & all(result$x[1] >=subset_data1$Ci)){
        
        judge <- TRUE
        
        predicted_values1 <- predict(model1)
        predicted_values2 <- predict(model2)
        
        rmse1 <- sqrt(mean((subset_data1$A - predicted_values1)^2))
        rmse2 <- sqrt(mean((subset_data2$A - predicted_values2)^2))
        
        # T_RMSEを計算
        T_RMSE <- rmse1 + rmse2
      } else {
        judge <- FALSE
        T_RMSE <- NA
      }
      
      #結果をデータフレームに入れていく。（.はパイプの中の前のやつを挿す）
      Judge_df <-
        data.frame(kind = i, Rubisco = j, RuBP = nrow(ACi_data)-j, 判定 = judge,　RMSE = T_RMSE) %>% 
        bind_rows(Judge_df, .)
      
      #結果出力
      print(Judge_df)
      
    }
    
    #データをRubisco回帰とRuBP回帰をするデータで分ける。
    Judging <- Judge_df %>% arrange(RMSE)
    
    if(Judging$判定[[1]] == "TRUE"){##########################################################################################################
      
      pre_data <- ACi_data %>% select(A, Ci)
      Rubisco_data <- pre_data[1:Judging[1,2],]
      RuBP_data <- pre_data[nrow(pre_data)-as.numeric(Judging[1,3])+1:nrow(pre_data),]
      RuBP_data <- na.omit(RuBP_data)
      Ci_range <- seq(0, 1500, by = 1)
      
      
      ##Rubisco回帰
      Rubisco_model <- nlsLM(A ~ Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd,
                             data = Rubisco_data,
                             start=c(Vcmax=20, Rd=0.2),
                             lower = c(Vcmax=6, Rd=0),
                             control = nls.control(maxiter = 100))
      Vcmax <- coef(Rubisco_model)[1]
      Rd <- coef(Rubisco_model)[2]
      Vcmax_SE <- summary(Rubisco_model)$coefficients[1, "Std. Error"]
      
      model_function_1 <- function(Ci, Vcmax, Rd){
        result_1 <- Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd
        return(result_1)
      }
      
      y_1 <- model_function_1(Ci_range, Vcmax, Rd)
      data_1 <- data.frame(Ci = Ci_range, A = y_1)
      
      ##RuBP回帰
      RuBP_model <- nlsLM(A ~ Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd, data = RuBP_data, start=c(Jmax=10))
      Jmax <- coef(RuBP_model)[1]
      
      
      model_function_2 <- function(Ci, Jmax){
        result_2 <- Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd
        return(result_2)
      }
      
      y_2 <- model_function_2(Ci_range, Jmax)
      data_2 <- data.frame(Ci = Ci_range, A = y_2)
      
      ##Ci_transitionを出す
      difference_function <- function(Ci, Vcmax, Rd, Jmax) {
        result_1 <- Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd
        result_2 <- Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd
        return(result_1 - result_2)
      }
      
      
      # 二つの関数の差が0となるCiの範囲を求める
      intersection <- uniroot(difference_function, interval = c(100, 1500), Vcmax = Vcmax, Rd = Rd, Jmax = Jmax)
      # 交点のCi
      intersection_Ci <- intersection$root
      # 交点のA
      intersection_A <- model_function_1(intersection_Ci, Vcmax, Rd)
      
      Ci_transition <- data.frame(Ci = intersection_Ci, A = intersection_A)
      
      # under limitの線
      filtered_data_1 <- data_1 %>% 
        filter(Ci <= Ci_transition$Ci)
      filtered_data_2 <- data_2 %>% 
        filter(Ci >= Ci_transition$Ci)
      
      # over limitの線
      out_data_1 <- data_1 %>% 
        filter(Ci > Ci_transition$Ci)
      out_data_2 <- data_2 %>% 
        filter(Ci < Ci_transition$Ci)
      
      ##plot
      plot <-
        ggplot()+
        geom_line(data = filtered_data_1, aes(x=Ci, y=A, color = "Rubisco limitation"), linewidth=1.5, show.legend = F)+
        geom_line(data = filtered_data_2, aes(x=Ci, y=A, color = "RuBP limitation"), linewidth=1.5, show.legend = F)+
        geom_line(data = out_data_1, aes(x=Ci, y=A), color = "#ff9900", linewidth=1.5, show.legend = F, alpha = 0.2)+
        geom_line(data = out_data_2, aes(x=Ci, y=A), color = "#339900", linewidth=1.5, show.legend = F, alpha = 0.2)+
        scale_color_manual(values = c("Rubisco limitation" = "#ff9900",
                                      "RuBP limitation" = "#339900"))+
        geom_point(data = ACi_data, aes(x=Ci, y=A))+
        geom_point(data = Ci_transition, aes(x=Ci, y=A), size=3.0, fill = "#ffd700", shape = 21)+
        labs(title = plant_name,
             x = expression(paste(italic(C)[i], " (", mu*mol, " ", {mol}^-1, ")")),
             y = expression(paste(italic(A), " (", mu*mol, " ", {{m}^-2}, " ",{s}^-1, ")")))+
        scale_x_continuous(breaks = seq(0, 1500, by = 300))+
        scale_y_continuous(breaks = seq(-5, 30, by =5))+
        xlim(0,1500)+
        ylim(-5,30)+
        annotate("text",
                 x = intersection_Ci,
                 y = intersection_A,
                 label = paste("(", round(intersection_Ci, 2), ",", round(intersection_A, 2), ")"),
                 vjust = -2,
                 hjust = 0.07,
                 color = "black")+
        theme(legend.position = "bottom", legend.justification = "center") +
        theme_classic(base_size = 15)
      
      
      #最終的にまとめる
      summarise_df <-
        data.frame(kind = plant_name, Rubisco = Judging[1,2], RuBP = Judging[1,3], RMSE = Judging[1,5],　Vcmax25 = Vcmax/Vc(t), Jmax25 = Jmax/J(t), Rd25 = Rd/Rd_25(t), VcmaxSE = Vcmax_SE) %>%
        bind_rows(summarise_df, .)
      
      plot_list[[i]] <- plot
    }
    else{ #######################################################################################################################
      
      pre_data <- ACi_data %>% select(A, Ci)
      Rubisco_data <- pre_data[1:5,]
      RuBP_data <- pre_data[5+1:nrow(pre_data),]
      RuBP_data <- na.omit(RuBP_data)
      Ci_range <- seq(0, 1500, by = 1)
      
      
      ##Rubisco回帰
      Rubisco_model <- nlsLM(A ~ Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd,
                             data = Rubisco_data,
                             start=c(Vcmax=20, Rd=0.2),
                             lower = c(Vcmax=6, Rd=0),
                             control = nls.control(maxiter = 100))
      Vcmax <- coef(Rubisco_model)[1]
      Rd <- coef(Rubisco_model)[2]
      Vcmax_SE <- summary(Rubisco_model)$coefficients[1, "Std. Error"]
      
      model_function_1 <- function(Ci, Vcmax, Rd){
        result_1 <- Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd
        return(result_1)
      }
      
      y_1 <- model_function_1(Ci_range, Vcmax, Rd)
      data_1 <- data.frame(Ci = Ci_range, A = y_1)
      
      ##RuBP回帰
      RuBP_model <- nlsLM(A ~ Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd, data = RuBP_data, start=c(Jmax=10))
      Jmax <- coef(RuBP_model)[1]
      
      
      model_function_2 <- function(Ci, Jmax){
        result_2 <- Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd
        return(result_2)
      }
      
      y_2 <- model_function_2(Ci_range, Jmax)
      data_2 <- data.frame(Ci = Ci_range, A = y_2)
      
      ##Ci_transitionを出す
      difference_function <- function(Ci, Vcmax, Rd, Jmax) {
        result_1 <- Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd
        result_2 <- Jmax*(Ci-g(t))/(4*Ci+8*g(t))-Rd
        return(result_1 - result_2)
      }
      
      
      # 二つの関数の差が0となるCiの範囲を求める
      intersection <- uniroot(difference_function, interval = c(100,500), Vcmax = Vcmax, Rd = Rd, Jmax = Jmax)
      # 交点のCi
      intersection_Ci <- intersection$root
      # 交点のA
      intersection_A <- model_function_1(intersection_Ci, Vcmax, Rd)
      
      Ci_transition <- data.frame(Ci = intersection_Ci, A = intersection_A)
      
      # under limitの線
      filtered_data_1 <- data_1 %>% 
        filter(Ci <= Ci_transition$Ci)
      filtered_data_2 <- data_2 %>% 
        filter(Ci >= Ci_transition$Ci)
      
      # over limitの線
      out_data_1 <- data_1 %>% 
        filter(Ci > Ci_transition$Ci)
      out_data_2 <- data_2 %>% 
        filter(Ci < Ci_transition$Ci)
      
      ##plot
      plot <-
        ggplot()+
        geom_line(data = filtered_data_1, aes(x=Ci, y=A, color = "Rubisco limitation"), linewidth=1.5, show.legend = F)+
        geom_line(data = filtered_data_2, aes(x=Ci, y=A, color = "RuBP limitation"), linewidth=1.5, show.legend = F)+
        geom_line(data = out_data_1, aes(x=Ci, y=A), color = "#ff9900", linewidth=1.5, show.legend = F, alpha = 0.2)+
        geom_line(data = out_data_2, aes(x=Ci, y=A), color = "#339900", linewidth=1.5, show.legend = F, alpha = 0.2)+
        scale_color_manual(values = c("Rubisco limitation" = "#ff9900",
                                      "RuBP limitation" = "#339900"))+
        geom_point(data = ACi_data, aes(x=Ci, y=A))+
        geom_point(data = Ci_transition, aes(x=Ci, y=A), size=3.0, fill = "#ffd700", shape = 21)+
        labs(title = plant_name,
             x = expression(paste(italic(C)[i], " (", mu*mol, " ", {mol}^-1, ")")),
             y = expression(paste(italic(A), " (", mu*mol, " ", {{m}^-2}, " ",{s}^-1, ")")))+
        scale_x_continuous(breaks = seq(0, max(filtered_data_2$Ci+100), by = 300))+
        scale_y_continuous(breaks = seq(-5, max(filtered_data_2$A+5), by =5))+
        xlim(0,max(ACi_data$Ci+100))+
        ylim(-5,max(ACi_data$A+5))+
        annotate("text",
                 x = intersection_Ci,
                 y = intersection_A,
                 label = paste("(", round(intersection_Ci, 2), ",", round(intersection_A, 2), ")"),
                 vjust = -2,
                 hjust = 0.07,
                 color = "black")+
        theme(legend.position = "bottom", legend.justification = "center") +
        theme_classic(base_size = 15)
      
      
      #最終的にまとめる
      summarise_df <-
        data.frame(kind = plant_name, Rubisco = 5, RuBP = nrow(pre_data)-5, RMSE = NA,　Vcmax25 = Vcmax/Vc(t), Jmax25 = Jmax/J(t), Rd25 = Rd/Rd_25(t), VcmaxSE = Vcmax_SE) %>%
        bind_rows(summarise_df, .)
      
      plot_list[[i]] <- plot
      
    }
    
    
  }
  
  
  
  
}



# Broad leaf --------------------------------------------------------------



# Graph -------------------------------------------------------------------



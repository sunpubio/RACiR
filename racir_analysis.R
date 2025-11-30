

# Library setting ---------------------------------------------------------
library(racir)
library(tidyverse)
library(ggplot2)
library(minpack.lm)
library(nleqslv)
library(Li6800fixer)

# Folder setting ----------------------------------------------------------
getwd()

broad_dir <- "data/broad"
conifer_dir <- "data/conifer"

list_broadleaf <- list.dirs(broad_dir, full.names = T, recursive = F)
list_conifer <- list.dirs(conifer_dir, full.names = T, recursive = F)



# df list setting --------------------------------------------------------------
Judge_df <- tibble()
summarise_df <- tibble()
plot_list <- list()

# parameter setting -------------------------------------------------------
temp.res <- function(p25, C, dH, t) exp(C-dH*1000/((t+273)*8.314))
g <- function(t) temp.res(42.75, 19.02, 37.83, t)
Kc <- function(t) temp.res(404.9, 38.05, 79.43, t)
Ko<- function(t) temp.res(278.4, 30.30, 36.38, t)
Vc <- function(t) temp.res(1, 26.35, 65.33, t)
Vo <- function(t) temp.res(1, 22.98, 60.11, t)
J <- function(t) temp.res(1, 17.71, 43.9, t)
Rd_25 <- function(t) temp.res(1, 18.72, 46.49, t)
K <- function(t) Kc(t)/(1+200*1000/Ko(t))
O2 <- 210


# Conifer -----------------------------------------------------------------

for (i in list_conifer) {
  
  print(i)
  
  all_files <- list.files(i, full.names = T, recursive = F) %>% .[!grepl("\\.xlsx?$", .)]
  aci_files <- all_files[grepl("\\_ACi", all_files)]
  area_files <- all_files[grepl("\\.csv", all_files)]
  
  library(dplyr)
  
  area_file <- read.csv(area_files, header = TRUE) %>% 
    mutate(
      # ① 日付以降があれば削る（なければそのまま）
      leaf = sub("_[0-9]{8}_.*$", "", File),
      # ② 拡張子（.jpg, .png など）を削る
      leaf = sub("\\.[^.]+$", "", leaf),
      # ③ mm^2 → m^2
      area_m2 = Area_mm2 / 1e6
    )
  
  
  for (m in aci_files) {
    
    print(m)
    
    ACi_data <- read_6800(m)
    plant_name <- str_extract(m, "[^/]+(?=_ACi)")
    ACi_data$PlantID <- plant_name
    
    area <- area_file %>% 
      filter(leaf == plant_name) %>% 
      pull(area_m2)
    
    fixed_ACi_data <- fixarea_6800(ACi_data, area)
    
    
    # 葉温決定
    t <- mean(ACi_data$Tleaf)
    
    for (j in 2:(nrow(fixed_ACi_data)-1)) {　　#点2つじゃ回帰出来なかったため、3つから行う。
      subset_data1 <- fixed_ACi_data[1:j, ] # 最初からi番目までのデータを取得
      subset_data2 <- fixed_ACi_data[(j + 1):nrow(fixed_ACi_data), ] # i+1番目から最後までのデータを取得
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
        data.frame(kind = i, Rubisco = j, RuBP = nrow(fixed_ACi_data)-j, 判定 = judge,　RMSE = T_RMSE) %>% 
        bind_rows(Judge_df, .)
      
      #結果出力
      print(Judge_df)
      
    }
    
    #データをRubisco回帰とRuBP回帰をするデータで分ける。
    Judging <- Judge_df %>% arrange(RMSE)
    
    if(Judging$判定[[1]] == "TRUE"){##########################################################################################################
      
      pre_data <- fixed_ACi_data %>% select(A, Ci)
      Rubisco_data <- pre_data[1:as.numeric(Judging[1,2]),]
      RuBP_data <- pre_data[nrow(pre_data)-as.numeric(as.numeric(Judging[1,3]))+1:nrow(pre_data),]
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
        geom_point(data = fixed_ACi_data, aes(x=Ci, y=A))+
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
      
      pre_data <- fixed_ACi_data %>% select(A, Ci)
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
        geom_point(data = fixed_ACi_data, aes(x=Ci, y=A))+
        geom_point(data = Ci_transition, aes(x=Ci, y=A), size=3.0, fill = "#ffd700", shape = 21)+
        labs(title = plant_name,
             x = expression(paste(italic(C)[i], " (", mu*mol, " ", {mol}^-1, ")")),
             y = expression(paste(italic(A), " (", mu*mol, " ", {{m}^-2}, " ",{s}^-1, ")")))+
        scale_x_continuous(breaks = seq(0, max(filtered_data_2$Ci+100), by = 300))+
        scale_y_continuous(breaks = seq(-5, max(filtered_data_2$A+5), by =5))+
        xlim(0,max(fixed_ACi_data$Ci+100))+
        ylim(-5,max(fixed_ACi_data$A+5))+
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

print(summarise_df)

# RACiR conifer -----------------------------------------------------------

# file setting
file_path <- "C:/Users/sabo/OneDrive - The University of Tokyo/研究/光合成能力推定手法/RACiR/test/data/conifer"

patterns <- c("50", "100", "200", "300", "400")
plot_racir_list <- list()
racir_df <- data.frame()
cal_time_plot <- list()


for (i in list_conifer) {
  
  
# どの葉なのか
  print(i)
  
  all_files <- list.files(i, full.names = T, recursive = F) %>% .[!grepl("\\.xlsx?$", .)]
  area_files <- all_files[grepl("\\.csv", all_files)]
  
　# 50, 100, 200, 300, 400で繰り返す
  for (pattern in patterns) {
    
    # patternに合うファイルを選ぶ
    matching_files <- all_files[grepl(pattern, all_files)]
    
    # calibrationとdataのファイルを分ける
    cal_der <- matching_files[grepl("\\_cal", matching_files)]
    data_der <- matching_files[-grepl("\\_cal", matching_files)]

    # calのセッティング
    caldata <- read_6800(cal_der)
    names(caldata) <- make.unique(names(caldata))
    n <- nrow(caldata)
    
    # 面積を修正
    area <- area_file %>% 
      filter(leaf == plant_name) %>% 
      pull(area_m2)
    
    caldata <- fixarea_6800(caldata, area)
    
    
    caldata$delta <- NA
    
    if (n > 1) {  # 行数が2以上であることを確認
      caldata$delta[2:n] <- caldata$A[2:n] - caldata$A[1:(n - 1)]
    }
    
    # calibrationの回帰に使うデータを選ぶ
    filtered_data <- caldata[abs(caldata$delta) <= 0.01 & caldata$CO2_r <= 400 & caldata$CO2_r >=50, c("obs", "CO2_r", "A", "delta")]

    
    ## ここから作図のための処理-------------------------------------------------------------
    　　cal <- caldata %>% 
      　　select(obs, delta, CO2_r, CO2_s, A, hhmmss)
    
    　　cal$hhmmss <- hms(cal$hhmmss)
    
    # obs=1の時刻を基準（0秒）に設定
    　　cal$seconds <- as.numeric(cal$hhmmss - cal$hhmmss[cal$obs == 1])
    
    #空のチャンバーでCO2を下げているデータ
    　　{ggplot()+
        geom_point(data =cal, aes(x=CO2_r, y=A))}
    
    #空のチャンバーでCO2_sとCO2_rで時間差があるという図を見ぜる
    {　　ggplot()+
    　   　 geom_point(data = cal, aes(x=seconds, y=CO2_r),
                   colour = "darkgreen")+
       　 geom_point(data = cal, aes(x=seconds, y=CO2_s),
                   colour= "lightgreen")+
        　labs(title = "Lags between the reference [CO2] and sample [CO2]",
             x= paste("time", " (s)"),
             y= expression(paste("[", italic(CO)[2], "]",  "  (", mu*mol, " ", {mol}^-1, ")")))+
       　 scale_color_manual(values = c("CO2_r" = "darkgreen", "CO2_s" = "lightgreen"))+
       　 theme_bw(base_size =16)}
    
   　　 df_long <- tidyr::pivot_longer(cal, cols = c(CO2_r, CO2_s), names_to = "type", values_to = "value")
    
    
    # ggplotでプロット
   　　 {cal_plot <- 
       　　 ggplot(df_long, aes(x = seconds, y = value, colour = type)) +
       　　 geom_point() +
       　　 labs(title = bquote(paste("Lags between the reference ", "[", italic(CO)[2],"]",  " and sample", "[", italic(CO)[2], "]", "  RR:", .(pattern))),
             x = "time (s)",
             y = expression(paste("[", italic(CO)[2], "]",  "  (", mu*mol, " ", {mol}^-1, ")")),
             colour = NULL) +
       　　 scale_color_manual(labels = c(expression(paste(italic(CO)[2], "_r")), expression(paste(italic(CO)[2], "_s"))),
                           values = c("darkgreen", "lightgreen")) +
       　　 theme_light(base_size = 16)+
       　　 theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 2),
              legend.position = c(0.85, 0.85))}
    
    
    #Ai-Ai-1をしてδにしたデータ
    　　{ggplot()+
       　　 geom_point(data = cal, aes(x=CO2_r, y=delta))+
       　　 geom_point(data = filtered_data, aes(x=CO2_r, y=delta, alpha=0.1), colour = "red", show.legend = F)+
       　　 ggtitle("filterd_delta")}
   
   　　 ggplot()+
      　　geom_point(data = cal, aes(x=CO2_r, y=A))+
     　　 geom_point(data = filtered_data, aes(x=CO2_r, y=A), colour = "red")
   
   ## ここまで作図----------------------------------------------------------------------------
    
   
   # キャリブレーション用のモデルを回す
    quadratic_model <- lm(A ~ poly(CO2_r, 1, raw = TRUE), data = filtered_data)
    
    
    # モデルが作成されたかどうかを確認
    if (nrow(filtered_data) <2) {
      message("エラーが発生しました。処理をスキップします。")
      
      # モデルが作れなかった時の処理
      for (m in seq_along(data_der)) {
        print(m)
        file <- data_der[m]
        data_data <- read_6800(file)
        names(data_data) <- make.unique(names(data_data))
        
        plant_name <- str_extract(
          file,
          paste0("[^/]+(?=_", pattern, "$)")  # ここで ) を足す
        )
        
        racir_df <- data.frame(Plantname = plant_name, RampRate = pattern, Method="RACiR" ,Vcmax25 = NA, SE = NA) %>% 
          bind_rows(racir_df,.)
        
        next}
     } else {
        
        # モデルが正常に作成された場合の処理
        for (m in seq_along(data_der)) {
          file <- data_der[m]
          data_data <- read_6800(file)
          names(data_data) <- make.unique(names(data_data))
          
          # 面積を修正
          area <- area_file %>% 
            filter(leaf == plant_name) %>% 
            pull(area_m2)
          
          data_data <- fixarea_6800(data_data, area)
          
          plant_name <- str_extract(
            file,
            paste0("[^/]+(?=_", pattern, "$)")  # ここで ) を足す
          )
          
          # calibrationデータで取ったCO2_rの範囲のデータを取得
          data_leaf <- data_data %>% 
            select(obs, A, Ci, CO2_r, Tleaf, E, gtc, Ca) %>% 
            filter(CO2_r >= min(filtered_data$CO2_r, na.rm = T) & CO2_r<=max(filtered_data$CO2_r, na.rm = T))   
          #filter(CO2_r <=380 & CO2_r>=90) %>% 
          #filter(obs >=10)
          
          
          # calの式で補正
          predicted_A <- predict(quadratic_model, newdata = data.frame(CO2_r = data_leaf$CO2_r))
          data_leaf$Corr <- predicted_A
          data_leaf$Aleaf <- data_leaf$A-data_leaf$Corr
          data_leaf$Ci <- ((data_leaf$gtc-(data_leaf$E/2))*data_leaf$Ca-data_leaf$Aleaf)/(data_leaf$gtc+data_leaf$E/2)
          data_leaf$A <- data_leaf$Aleaf
          
          t <- mean(data_leaf$Tleaf)
          
          area <- area_file %>% 
            filter(leaf == plant_name) %>% 
            pull(area_m2)
          
          fixed_ACi_data <- fixarea_6800(data_leaf, area)
          
          
          # Vcmax推定
          model1 <- nlsLM(A ~ Vcmax*(Ci-g(t))/(Ci+Kc(t)*(1+O2/Ko(t)))-Rd, data = data_leaf, start=c(Vcmax=50, Rd=0.8), control = nls.control(maxiter = 100))
          Vcmax <- coef(model1)[1]
          Rd <- coef(model1)[2]
          
          
          # デーらフレームに結果を保存
          racir_df <- data.frame(Plantname = plant_name, RampRate = pattern, Method="RACiR" ,Vcmax25 = Vcmax/Vc(t), SE = NA) %>% 
            bind_rows(racir_df,.)
          
          
          # グラフ保存
          racir_plot<-
            plot_list[[m]]+
            geom_point(data = data_leaf, aes(x=Ci, y=A), colour="red", alpha=0.6)+
            ggtitle(paste(plant_name, "  ACi  &  RACiR", pattern))
          
        
        } 
      }
    
    }
  }
  
  
  
  
  
  
  
  
  











# Broad leaf --------------------------------------------------------------



# Graph -------------------------------------------------------------------



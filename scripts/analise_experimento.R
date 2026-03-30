# ============================================================
# ANALISE ESTATISTICA - ahmdados.csv
# Estrutura esperada do projeto:
# analise-adubacao-nitrogenio/
# ├─ data/
# │  └─ ahmdados.csv
# ├─ scripts/
# │  └─ analise_experimento.R
# ├─ results/
# ├─ docs/
# │  └─ figuras/
# ============================================================

# ------------------------------
# 1) Pacotes
# ------------------------------
packs <- c("ggplot2", "dplyr", "emmeans", "broom", "car")
instalar <- packs[!sapply(packs, requireNamespace, quietly = TRUE)]
if (length(instalar) > 0) install.packages(instalar, dependencies = TRUE)
invisible(lapply(packs, library, character.only = TRUE))

# ------------------------------
# 2) Diretorios do projeto
# ------------------------------
# Assumindo que o script sera executado a partir da pasta scripts/
dir_scripts <- getwd()
dir_projeto <- normalizePath(file.path(dir_scripts, ".."), winslash = "/", mustWork = TRUE)

dir_data    <- file.path(dir_projeto, "data")
dir_results <- file.path(dir_projeto, "results")
dir_docs    <- file.path(dir_projeto, "docs")
dir_figuras <- file.path(dir_docs, "figuras")

if (!dir.exists(dir_results)) dir.create(dir_results, recursive = TRUE)
if (!dir.exists(dir_docs)) dir.create(dir_docs, recursive = TRUE)
if (!dir.exists(dir_figuras)) dir.create(dir_figuras, recursive = TRUE)

arquivo <- file.path(dir_data, "ahmdados.csv")

if (!file.exists(arquivo)) {
  stop(paste("Arquivo nao encontrado em:", arquivo))
}

# ------------------------------
# 3) Leitura segura
# ------------------------------
dados_brutos <- read.csv(
  arquivo,
  header = TRUE,
  sep = ",",
  dec = ".",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------
# 4) Inspecao das colunas
# ------------------------------
resumo_colunas <- data.frame(
  coluna = names(dados_brutos),
  classe = sapply(dados_brutos, function(x) class(x)[1]),
  n_unicos = sapply(dados_brutos, function(x) length(unique(x)))
)

print("Resumo das colunas:")
print(resumo_colunas)

cols_necessarias <- c("ADU", "LIN", "REP", "PROD")
faltantes <- setdiff(cols_necessarias, names(dados_brutos))
if (length(faltantes) > 0) {
  stop(paste("Colunas ausentes no arquivo:", paste(faltantes, collapse = ", ")))
}

# ------------------------------
# 5) Padronizacao de tipos
# ------------------------------
dados <- dados_brutos

dados$ADU_num <- suppressWarnings(as.numeric(dados$ADU))
dados$ADU     <- factor(dados$ADU)
dados$LIN     <- factor(dados$LIN)
dados$REP     <- factor(dados$REP)
dados$PROD    <- suppressWarnings(as.numeric(dados$PROD))

dados$LIN_nome <- factor(
  dados$LIN,
  levels = levels(dados$LIN),
  labels = c("Simples", "Duplas")[seq_along(levels(dados$LIN))]
)

# ------------------------------
# 6) Delineamento identificado
# ------------------------------
delineamento_txt <- paste(
  "DBC fatorial 4 x 2 com 3 blocos, considerando ADU como fator quantitativo,",
  "LIN como fator qualitativo, REP como bloco e PROD como variavel resposta."
)

# ------------------------------
# 7) ANOVA fatorial em DBC
# ------------------------------
mod_fatorial <- aov(PROD ~ REP + ADU * LIN, data = dados)
anova_fatorial <- anova(mod_fatorial)

# ------------------------------
# 8) Residuos
# ------------------------------
res <- residuals(mod_fatorial)

shapiro_res <- shapiro.test(res)
levene_res <- car::leveneTest(PROD ~ interaction(ADU, LIN), data = dados)

# ------------------------------
# 9) Medias por combinacao fatorial
# ------------------------------
medias <- dados %>%
  group_by(ADU, LIN, LIN_nome) %>%
  summarise(
    n = n(),
    media = mean(PROD, na.rm = TRUE),
    dp = sd(PROD, na.rm = TRUE),
    ep = dp / sqrt(n),
    .groups = "drop"
  )

# ------------------------------
# 10) Desdobramento de LIN dentro de ADU
# ------------------------------
emm_lin_dentro_adu <- emmeans(mod_fatorial, ~ LIN | ADU)
compar_lin_dentro_adu <- contrast(
  emm_lin_dentro_adu,
  method = "pairwise",
  adjust = "tukey"
) %>%
  broom::tidy()

# ------------------------------
# 11) Regressao de ADU dentro de cada LIN (dados brutos)
# ------------------------------
lin_niveis <- levels(dados$LIN)

mod_lin1 <- lm(PROD ~ REP + ADU_num + I(ADU_num^2),
               data = subset(dados, LIN == lin_niveis[1]))

mod_lin2 <- lm(PROD ~ REP + ADU_num + I(ADU_num^2),
               data = subset(dados, LIN == lin_niveis[2]))

anova_lin1 <- anova(mod_lin1)
anova_lin2 <- anova(mod_lin2)

coef_lin1 <- summary(mod_lin1)
coef_lin2 <- summary(mod_lin2)

# ------------------------------
# 12) Comparacao marginal de LIN
# ------------------------------
emm_lin <- emmeans(mod_fatorial, ~ LIN)
compar_lin <- contrast(emm_lin, method = "pairwise") %>%
  broom::tidy()

# ------------------------------
# 13) Tabelas para exportacao
# ------------------------------
tab_anova <- broom::tidy(mod_fatorial)
tab_reg_lin1 <- broom::tidy(mod_lin1)
tab_reg_lin2 <- broom::tidy(mod_lin2)

tab_shapiro <- data.frame(
  secao = "shapiro_residuos",
  teste = "Shapiro-Wilk",
  estatistica = unname(shapiro_res$statistic),
  p_valor = shapiro_res$p.value
)

tab_levene <- data.frame(
  secao = "levene",
  teste = "Levene",
  estatistica = levene_res$`F value`[1],
  p_valor = levene_res$`Pr(>F)`[1]
)

saida_tabelas <- bind_rows(
  tab_anova %>% mutate(secao = "anova_fatorial"),
  medias %>% mutate(secao = "medias_tratamentos"),
  compar_lin_dentro_adu %>% mutate(secao = "comparacao_LIN_dentro_ADU"),
  compar_lin %>% mutate(secao = "comparacao_marginal_LIN"),
  tab_reg_lin1 %>% mutate(secao = "regressao_LIN1_bruto"),
  tab_reg_lin2 %>% mutate(secao = "regressao_LIN2_bruto"),
  tab_shapiro,
  tab_levene
)

write.csv(
  saida_tabelas,
  file.path(dir_results, "tabelas_resultados.csv"),
  row.names = FALSE
)

# ------------------------------
# 14) Tema de grafico
# ------------------------------
theme_artigo <- theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey85"),
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_text(face = "bold")
  )

# ------------------------------
# 15) Grafico de barras com EP
# ------------------------------
grafico_barras <- ggplot(medias, aes(x = ADU, y = media, fill = LIN_nome)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = media - ep, ymax = media + ep),
    position = position_dodge(width = 0.8),
    width = 0.2,
    linewidth = 0.6
  ) +
  labs(
    title = "Produtividade por sistema de fileira dentro de cada nivel de adubacao",
    x = "Adubacao (kg/ha de Nitrogenio)",
    y = expression("Produtividade (kg ha"^{-1}*")"),
    fill = "Sistema de fileira"
  ) +
  theme_artigo

ggsave(
  file.path(dir_results, "grafico_barras_EP.png"),
  plot = grafico_barras,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(dir_figuras, "grafico_barras_EP.png"),
  plot = grafico_barras,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------
# 16) Grafico de interacao
# ------------------------------
grafico_interacao <- ggplot(
  medias,
  aes(x = as.numeric(as.character(ADU)), y = media, color = LIN_nome, group = LIN_nome)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = media - ep, ymax = media + ep), width = 3, linewidth = 0.7) +
  labs(
    title = "Interacao entre adubacao e sistema de fileira",
    x = "Adubacao (kg/ha de Nitrogenio)",
    y = "Produtividade (media ± EP)",
    color = "Sistema de fileira"
  ) +
  theme_artigo

ggsave(
  file.path(dir_results, "interacao_ADU_LIN.png"),
  plot = grafico_interacao,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(dir_figuras, "interacao_ADU_LIN.png"),
  plot = grafico_interacao,
  width = 8,
  height = 6,
  dpi = 300
)

# ------------------------------
# 17) Regressao com dados brutos por LIN
# ------------------------------
pred1 <- data.frame(
  ADU_num = seq(min(dados$ADU_num), max(dados$ADU_num), length.out = 200),
  REP = levels(dados$REP)[1]
)
pred1$fit <- predict(mod_lin1, newdata = pred1)

dados_lin1 <- subset(dados, LIN == lin_niveis[1])
grafico_reg1 <- ggplot(dados_lin1, aes(x = ADU_num, y = PROD)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(data = pred1, aes(x = ADU_num, y = fit), linewidth = 1) +
  labs(
    title = "Regressao de PROD em ADU - LIN 1",
    x = "Adubacao (kg/ha de Nitrogenio)",
    y = "Produtividade"
  ) +
  theme_artigo

ggsave(
  file.path(dir_results, "regressao_LIN1_bruto.png"),
  plot = grafico_reg1,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(dir_figuras, "regressao_LIN1_bruto.png"),
  plot = grafico_reg1,
  width = 8,
  height = 6,
  dpi = 300
)

pred2 <- data.frame(
  ADU_num = seq(min(dados$ADU_num), max(dados$ADU_num), length.out = 200),
  REP = levels(dados$REP)[1]
)
pred2$fit <- predict(mod_lin2, newdata = pred2)

dados_lin2 <- subset(dados, LIN == lin_niveis[2])
grafico_reg2 <- ggplot(dados_lin2, aes(x = ADU_num, y = PROD)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_line(data = pred2, aes(x = ADU_num, y = fit), linewidth = 1) +
  labs(
    title = "Regressao de PROD em ADU - LIN 2",
    x = "Adubacao (kg/ha de Nitrogenio)",
    y = "Produtividade"
  ) +
  theme_artigo

ggsave(
  file.path(dir_results, "regressao_LIN2_bruto.png"),
  plot = grafico_reg2,
  width = 8,
  height = 6,
  dpi = 300
)

ggsave(
  file.path(dir_figuras, "regressao_LIN2_bruto.png"),
  plot = grafico_reg2,
  width = 8,
  height = 6,
  dpi = 300
)

# ============================================================
# 18) REGRESSAO LINEAR E QUADRATICA COM MEDIAS DOS TRATAMENTOS
# ============================================================

dados_medias <- medias %>%
  transmute(
    ADU = as.numeric(as.character(ADU)),
    LIN = as.character(LIN_nome),
    PROD = media
  )

calc_r2 <- function(modelo, y) {
  y_aj <- fitted(modelo)
  ss_res <- sum((y - y_aj)^2)
  ss_tot <- sum((y - mean(y))^2)
  1 - (ss_res / ss_tot)
}

equacao_linear <- function(mod) {
  b <- coef(mod)
  paste0(
    "y = ", round(b[1], 2),
    ifelse(b[2] >= 0, " + ", " - "),
    abs(round(b[2], 4)), "x"
  )
}

equacao_quadratica <- function(mod) {
  b <- coef(mod)
  paste0(
    "y = ", round(b[1], 2),
    ifelse(b[2] >= 0, " + ", " - "),
    abs(round(b[2], 4)), "x",
    ifelse(b[3] >= 0, " + ", " - "),
    abs(round(b[3], 6)), "x²"
  )
}

resumo_saida <- c()
resultados <- list()
predicoes <- data.frame()

for (lin in levels(factor(dados_medias$LIN, levels = c("Simples", "Duplas")))) {
  
  sub <- subset(dados_medias, LIN == lin)
  
  mod_lin <- lm(PROD ~ ADU, data = sub)
  mod_quad <- lm(PROD ~ ADU + I(ADU^2), data = sub)
  
  r2_lin  <- calc_r2(mod_lin, sub$PROD)
  r2_quad <- calc_r2(mod_quad, sub$PROD)
  
  coefs <- coef(mod_quad)
  a1 <- coefs[2]
  a2 <- coefs[3]
  
  if (!is.na(a2) && a2 != 0) {
    met_x <- -a1 / (2 * a2)
    met_y <- as.numeric(predict(mod_quad, newdata = data.frame(ADU = met_x)))
  } else {
    met_x <- NA
    met_y <- NA
  }
  
  novo <- data.frame(ADU = seq(min(sub$ADU), max(sub$ADU), length.out = 200))
  pred <- predict(mod_quad, newdata = novo, se.fit = TRUE)
  t_crit <- qt(0.975, df = df.residual(mod_quad))
  
  novo$fit <- pred$fit
  novo$se  <- pred$se.fit
  novo$lwr <- novo$fit - t_crit * novo$se
  novo$upr <- novo$fit + t_crit * novo$se
  novo$LIN <- lin
  
  predicoes <- rbind(predicoes, novo)
  
  resultados[[lin]] <- list(
    r2_lin = r2_lin,
    r2_quad = r2_quad,
    eq_lin = equacao_linear(mod_lin),
    eq_quad = equacao_quadratica(mod_quad),
    met_x = met_x,
    met_y = met_y,
    anova_lin = capture.output(anova(mod_lin)),
    anova_quad = capture.output(anova(mod_quad)),
    resumo_quad = capture.output(summary(mod_quad))
  )
  
  resumo_saida <- c(
    resumo_saida,
    "=========================================================",
    paste("SISTEMA DE FILEIRA:", lin),
    "=========================================================",
    "",
    "Modelo linear:",
    paste("Equacao:", resultados[[lin]]$eq_lin),
    paste("R2 =", round(r2_lin, 4)),
    "",
    "Modelo quadratico:",
    paste("Equacao:", resultados[[lin]]$eq_quad),
    paste("R2 =", round(r2_quad, 4)),
    paste("MET (kg/ha de N) =", round(met_x, 2)),
    paste("Produtividade estimada no MET (kg/ha) =", round(met_y, 2)),
    "",
    "ANOVA do modelo linear:",
    resultados[[lin]]$anova_lin,
    "",
    "ANOVA do modelo quadratico:",
    resultados[[lin]]$anova_quad,
    "",
    "Resumo do modelo quadratico:",
    resultados[[lin]]$resumo_quad,
    ""
  )
}

writeLines(resumo_saida, file.path(dir_results, "relatorio_regressao_medias.txt"))

tabela_modelos <- data.frame(
  LIN = c("Simples", "Duplas"),
  Equacao_linear = c(resultados$Simples$eq_lin, resultados$Duplas$eq_lin),
  R2_linear = c(resultados$Simples$r2_lin, resultados$Duplas$r2_lin),
  Equacao_quadratica = c(resultados$Simples$eq_quad, resultados$Duplas$eq_quad),
  R2_quadratica = c(resultados$Simples$r2_quad, resultados$Duplas$r2_quad),
  MET_kg_ha_N = c(resultados$Simples$met_x, resultados$Duplas$met_x),
  PROD_no_MET_kg_ha = c(resultados$Simples$met_y, resultados$Duplas$met_y)
)

write.csv(
  tabela_modelos,
  file = file.path(dir_results, "tabela_regressao_medias.csv"),
  row.names = FALSE
)

anotacoes <- data.frame(
  LIN = c("Simples", "Duplas"),
  x = c(32, 32),
  y = c(
    max(predicoes$upr[predicoes$LIN == "Simples"]) * 0.98,
    max(predicoes$upr[predicoes$LIN == "Duplas"]) * 0.88
  ),
  label = c(
    paste0(
      "Simples\n",
      resultados$Simples$eq_quad, "\n",
      "R² = ", round(resultados$Simples$r2_quad, 3), "\n",
      "MET = ", round(resultados$Simples$met_x, 2), " kg/ha"
    ),
    paste0(
      "Duplas\n",
      resultados$Duplas$eq_quad, "\n",
      "R² = ", round(resultados$Duplas$r2_quad, 3), "\n",
      "MET = ", round(resultados$Duplas$met_x, 2), " kg/ha"
    )
  )
)

grafico_met <- ggplot() +
  geom_ribbon(
    data = predicoes,
    aes(x = ADU, ymin = lwr, ymax = upr, fill = LIN),
    alpha = 0.18
  ) +
  geom_line(
    data = predicoes,
    aes(x = ADU, y = fit, color = LIN),
    linewidth = 1.1
  ) +
  geom_point(
    data = dados_medias,
    aes(x = ADU, y = PROD, color = LIN),
    size = 3
  ) +
  geom_point(
    data = data.frame(
      ADU = c(resultados$Simples$met_x, resultados$Duplas$met_x),
      PROD = c(resultados$Simples$met_y, resultados$Duplas$met_y),
      LIN = c("Simples", "Duplas")
    ),
    aes(x = ADU, y = PROD, color = LIN),
    shape = 4, size = 4, stroke = 1.2
  ) +
  geom_text(
    data = anotacoes,
    aes(x = x, y = y, label = label, color = LIN),
    hjust = 0, vjust = 1, size = 3.6, show.legend = FALSE
  ) +
  labs(
    title = "Regressao quadratica da produtividade em funcao da adubacao nitrogenada",
    subtitle = "Ajuste com medias dos tratamentos para fileiras simples e duplas",
    x = "Adubacao (kg/ha de Nitrogenio)",
    y = expression("Produtividade (kg ha"^{-1}*")"),
    color = "Sistema de fileira",
    fill = "Sistema de fileira"
  ) +
  theme_artigo

ggsave(
  filename = file.path(dir_results, "regressao_quadratica_IC95_MET.tiff"),
  plot = grafico_met,
  width = 9,
  height = 6.5,
  dpi = 600,
  compression = "lzw"
)

# Análise Descritiva dos dados ####
# Leitura dos dados
rm(list = ls(all = TRUE))
dados <- read.csv2("lagartos.txt", sep = "", header = TRUE)
head(dados)
dados[,c("período", "comp.madeira", "larg.madeira", "ocupação")]<-lapply(dados[,c("período", "comp.madeira", "larg.madeira", "ocupação")], factor)
str(dados)
attach(dados)

#> Boxplot ####
library(png)
par(mfrow=c(1,2))
boxplot(grahani, col = "#cf3266", main = "Espécie Grahani", pch = 19, font = "serif")
boxplot(opalinus, col = "#cf3266", main = "Espécie Opalinus", pch = 19, font = "serif")


sort(dados$grahani)
sort(dados$opalinus)
dados[13,]
dados[16,]
dados[19,]

#> Medidas de Resumo ####
summary(dados[,1:2])

# Desvio padrão
dp.g <- sd(grahani); dp.g 
dp.o <- sd(opalinus); dp.o
sum(grahani)
sum(opalinus)

# Coeficiente de variação
cv.g <-  (dp.g/mean(grahani))*100; cv.g
cv.o <-  (dp.o/mean(opalinus))*100; cv.o

# resumo das variáveis categóricas
n <- length(grahani) # total amostral
names(dados)

# frequencia e proporçao para o período
f.per <- table(período); f.per
p.per <- rep((8/n)*100, 2); p.per
p.per2 <- rep((7/n)*100, 1); p.per2

# frequencia e proporçao para comp.madeira
f.cm <- table(comp.madeira); f.cm
p.cm <- (12/n)*100;p.cm
p.cm2 <- (11/n)*100;p.cm2

# frequencia e proporçao para larg.madeira
f.lm <- table(larg.madeira); f.lm
p.lm <- (12/n)*100; p.lm
p.lm2 <- (11/n)*100; p.lm2

# frequencia e proporçao para ocupação
f.o <- table(ocupação); f.o
p.o <- (12/n)*100; p.o
p.o2 <- (11/n)*100; p.o2

# associacao entre a resposta e as covariaveis categoricas
n.per <- c("manhã", "meio-dia", "tarde")
n.com <- c("curta", "comprida")
n.lar <- c("estreita", "larga")
n.ocu <- c("claro", "escuro")

par(mfrow=c(2,2))
boxplot(grahani~período, names=n.per, xlab="Período do dia", ylab="Espécie grahani", pch=19, col = "#cf3266")
boxplot(grahani~comp.madeira, names=n.com, xlab="Comprimento da madeira", ylab="Espécie grahani", pch=19, col = "#cf3266")
boxplot(grahani~larg.madeira, names=n.lar, xlab="Largura da madeira", ylab="Espécie grahani", pch=19, col = "#cf3266")
boxplot(grahani~ocupação, names=n.ocu, xlab="Local de ocupação", ylab="Espécie grahani", pch=19, col = "#cf3266")

# Ajuste do MLG ####
#> Caselas de referência ####
#| período: 1 - manhã
#| comp.madeira: 1 - curta
#| larg.madeira: 1 - estreita
#| ocupação: 1 - claro

# Variável resposta (binária)
Y <- cbind(grahani, opalinus)

# Modelo saturado sob H1 com interação de primeira ordem:
mod1 <- glm(Y ~ período*comp.madeira + período*larg.madeira + período*ocupação +
              comp.madeira*larg.madeira + comp.madeira*ocupação + larg.madeira * ocupação, 
            family = binomial(link = "logit"))
summary(mod1)
AIC(mod1)
BIC(mod1)
beta <- coef(mod1); beta # coeficientes estimados do MLG
odds <- exp(beta); odds  # chances estimadas

# Parâmetro de precisão esrimado phi
X <- model.matrix(mod1)
p <- ncol(X)
n <- nrow(X)
phi <- sum(residuals(mod1, type="pearson")^2)/(n-p); round(phi,4)

#> Seleção de modelo ####
library(MASS)
stepAIC(mod1, direction = "backward")
stepAIC(mod1)

# Modelo saturado sob H0:
mod0 <- glm(formula = Y ~ comp.madeira, family = binomial(link = "logit"))

AIC(mod0)
BIC(mod0)

#> Estimativas do Modelo ####
summary(mod0)

#> Qualidade do Ajuste ####
# Teste da Razão de verossimilhanças usando a funcao deviance do pacote stats
D1 <- deviance(mod1); D1
D0 <- deviance(mod0); D0
RV1 <- (D0 - D1); RV1
df <- df.residual(mod0) - df.residual(mod1); df
pvRV <- pchisq(RV1, df, lower.tail = FALSE); pvRV

# ou usando a funcao anova do pacote stats
RV2 <- anova(mod0, mod1,test='Chisq'); RV2
#| O teste não rejeitou o modelo restrito.

beta <- coef(mod0); beta # coeficientes estimados do MLG
odds <- exp(beta); odds  # chances estimadas

# Probabilidade de sucesso
mu <- exp(mod0$coefficients)/(1+exp(mod0$coefficients)); mu

# Parâmetro de precisão esrimado phi
X <- model.matrix(mod0)
p <- ncol(X)
n <- nrow(X)
phi <- sum(residuals(mod0, type="pearson")^2)/(n-p); round(phi,4)

### para os coeficientes
IC.beta <- confint(mod0, level=0.95); IC.beta

### para as chances
IC.odds <- exp(IC.beta); IC.odds

# Análise de Diagnóstico ####

X <- model.matrix(mod0)
p <- ncol(X)
n <- nrow(X)
eta <- mod0$fit
V <- fitted(mod0)
V <- diag(V)
w <- mod0$weights
W <- diag(w)
W1 <- sqrt(W)
H <- W1%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W1
h <- diag(H)
rd <- resid(mod0, type="deviance")
phi <- sum(residuals(mod0, type="pearson")^2)/(n-p); phi # precisao estimada
td <- rd*sqrt(phi/(1-h))
rp <- sqrt(phi)*resid(mod0, type="pearson")  # residuos de Pearson
ts <- rp/sqrt(1-h)                              # residuos de Pearson padronizados
ma <- max(td)
mi <- min(td)
LD <- h*(ts^2)/(1-h)                            # medidas de afastamento da verossimilhanca

#> Pontos de Alavanca ####
library(ggfortify)

data <- data.frame(index = 1:length(h), h = h)
alavanca <- hatvalues(mod0)

# Adicionar pontos de alavanca aos dados
dados$alavanca <- alavanca
#---------
ggplot(data, aes(x = index, y = h)) +
  geom_point(shape = 19) +
  geom_hline(yintercept = 2*p/n, linetype = "dashed", color = "#a30b3e", size = 1) +
  labs(x = "Índice", y = "Pontos de Alavanca", title = "") +
  ylim(0, 1) +
  theme_minimal() +
  geom_label(aes(label = ifelse(h > 2*(p/n), index, NA)), vjust = -1, hjust = 0.5)

#> Resíduos ####
## residuos do desvio - investigação de pontos aberrantes e heteroscedasticidade
datres <- data.frame(indx = 1:length(td), td = td)

ggplot(datres, aes(x = indx, y = td)) +
  geom_point(shape = 19) +
  geom_hline(yintercept = c(-2,2), linetype = "dashed", color = "#a30b3e", size = 1) +
  labs(x = "Índice", y = "Resíduos do Desvio", title = "") +
  ylim(-2.5, 2.5) +
  theme_minimal() 

#> Pontos Influentes ####
# Criando o gráfico com ggplot2
ggplot(data, aes(x = index, y = LD)) +
  geom_point(shape = 19) +
  ylim(0, 0.2) +
  geom_hline(yintercept = 2*p/n, linetype = "dashed", color = "#a30b3e", size = 1) +
  labs(x = "Índice", y = "Medida de Afastamento") +
  theme_minimal()

#> Normalidade dos Resíduos  ####
shapiro.test(td) # Resíduos do Desvio
#| Os resíduos apresentam normalidade.

shapiro.test(rp) # Resíduos de pearson
#| Os resíduos apresentam normalidade.

## funcao de autocorrelacao e densidade dos residuos do desvio
bacf <- acf(td, plot = FALSE)
bacfdf <- with(bacf, data.frame(lag, acf))

ggplot(data = bacfdf, aes(x = lag, y = acf)) +
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0)) +
  geom_hline(yintercept = 0.41, linetype = "dashed", color = '#a30b3e', size = 1) + 
  geom_hline(yintercept = -0.41, linetype = "dashed", color = '#a30b3e', size = 1)+
  theme_minimal()

ggplot(datres, aes(x = td)) +
  geom_density() +
  labs(x = "Resíduos do Desvio", y = "Densidade", title = "") +
  xlim(-0.8,0.8) + ylim(0,1.5) +
  theme_minimal()

#> Envelope ####
### Funcao Envelope da binomail com reposição
source("~/Github/MLGs/modelo-binomial_lagartos/Envelopes/envelope.binomialcr.R")

## Grafico de probabilidade normal com envelope
form <- Y ~comp.madeira
Fam <- binomial(link="logit")
ni <- round(mean(Y),0)
total <- matrix(ni, n, 1)

envelope.binomialcr(form = form, Fam = Fam, k = 100, alfa = 0.05, total = total)

# Análise Confirmatória ####
dados[13,]; dados[16,]; dados[19,]

#> Retirando observações: 
###  Sem a observacao 13 ###
ajuste.13 <- glm(Y ~ comp.madeira, subset = -c(13), family = binomial(link = "logit"))
summary(ajuste.13)

### Sem a observacao 16 ###
ajuste.16 <- glm(Y ~ comp.madeira, subset = -c(16), family = binomial(link = "logit"))
summary(ajuste.16)

### Sem a observação 19 ###
ajuste.19 <- glm(Y ~ comp.madeira, subset = -c(19), family = binomial(link = "logit"))
summary(ajuste.19)

### Sem as observações 13 e 16 ###
ajuste.1316 <- glm(Y ~ comp.madeira, subset = -c(13,16), family = binomial(link = "logit"))
summary(ajuste.1316)

### Sem as observações 13 e 19 ###
ajuste.1319 <- glm(Y ~ comp.madeira, subset = -c(13,19), family = binomial(link = "logit"))
summary(ajuste.1319)

### Sem as observações 16 e 19 ###
ajuste.1619 <- glm(Y ~ comp.madeira, subset = -c(16,19), family = binomial(link = "logit"))
summary(ajuste.1619)

### Sem as observações 13, 16 e 19 ###
ajuste.131619 <- glm(Y ~ comp.madeira, subset = -c(13,16,19), family = binomial(link = "logit"))
summary(ajuste.131619)

#########################################################################
# Variação Percentual das Estimativas do modelo sem as observações #
#########################################################################

VP.13 <- ((coef(ajuste.13)-coef(mod0))/coef(mod0))*100; VP.13
VP.16 <- ((coef(ajuste.16)-coef(mod0))/coef(mod0))*100; VP.16
VP.19 <- ((coef(ajuste.19)-coef(mod0))/coef(mod0))*100; VP.19
VP.1316 <- ((coef(ajuste.1316)-coef(mod0))/coef(mod0))*100; VP.1316
VP.1319 <- ((coef(ajuste.1319)-coef(mod0))/coef(mod0))*100; VP.1319
VP.1619 <- ((coef(ajuste.1619)-coef(mod0))/coef(mod0))*100; VP.1619
VP.131619 <- ((coef(ajuste.131619)-coef(mod0))/coef(mod0))*100; VP.131619

#------------------------------------------------------------------------------#
#### NOVO MODELO ####
#> Modelo quasi-binomial
mod00 <- glm(Y ~ comp.madeira, family = quasibinomial(link = "logit"))

hn <- hnp(mod00, resid.type = "pearson", conf = 0.95)

Gb <- hnp(hn$residuals,how.many.out = T,print.on = F, half = T, 
          resid.type = "pearson", conf = 0.95)
Gb1 <- hnp(hn$residuals,how.many.out = T,print.on = F, half = T, 
           resid.type = "pearson", conf = 0.90)
Gb2 <- hnp(hn$residuals,how.many.out = T,print.on = F, half = T, 
           resid.type = "pearson", conf = 0.99)

G1b <- with(Gb, data.frame(x, lower, upper, median, residuals))

#> Envelope para o modelo quasi-binomial
ggplot(data = G1b, aes(x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),fill="#cf3266",alpha=0.5)+
  geom_ribbon(aes(ymin=Gb1$lower, ymax=Gb1$upper),fill="#cf3266",alpha=0.5) +
  geom_ribbon(aes(ymin=Gb2$lower, ymax=Gb2$upper),fill="#cf3266",alpha=0.5) +
  scale_fill_gradient(low = "#a30b3e", high = "#fa87ad") + 
  geom_point(aes(y = residuals), show.legend = FALSE) + 
  geom_line(aes(y = median),lty=2) +
  xlab("Quantis Teóricos") + 
  ylab("Quantis Empíricos") + 
  theme(text=element_text(size=15,family="serif"))




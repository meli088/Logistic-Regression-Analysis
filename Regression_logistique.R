library("FactoMineR")
library("corrplot")

data <- read.table("mat.descriptors.dat", header = TRUE, sep = " ")
ncol_cont = ncol(data)-1
summary(data)

data_noNA <- na.omit(data)
summary(data_noNA)

variances <- sapply(data_noNA[-22], var)
variances

data_no12 <- data_noNA[,-21]
# Boxplot avant le nettoyage des données
boxplot(data_no12[1:20], las = 2,
        main = "Distribution des descripteurs avant nettoyage",
        xlab = "Descripteurs",
        ylab = "Valeurs")


xtrem_radius <- which(data_no12[,"radius"]>2700) # Identification de l'invididu
data_clean <- data_no12[-xtrem_radius,]

# Boxplot après le nettoyage des données
boxplot(data_clean[1:20], las = 2,
        main = "Distribution des descripteurs après nettoyage",
        xlab = "Descripteurs",
        ylab = "Valeurs")

# Effectuer l'analyse PCA sur les données nettoyées
pca <- PCA(data_clean[1:20], graph = T)

# Visualisation des résultats de PCA
plot(pca, choix = "ind", 
     title = "Analyse en composantes principales des molécules")


which(data_clean[1:20][, "a_nH"]>100) # 610  703 2849 2852 2854
which(data_clean[1:20][, "a_count"]>250) # 703

p.values <- sapply(data_clean[, -which(names(data_clean) == "drug")], function(x) {
  t.test(x ~ data_clean$drug)$p.value
})
p.values

seuil <- 0.05
variables_significatives <- names(p.values)[p.values < seuil]
variables_non_significatives <- names(p.values)[p.values >= seuil]

cat("Variables significatives (p < 0.05) :\n", variables_significatives, "\n\n")
cat("Variables non significatives (p >= 0.05) :\n", variables_non_significatives)

# Exclure nB et nBr car pas de différence significative entre dg et ndg.

data_clean_sans_non_signif <- data_clean[, -which(names(data_clean) %in% c("a_nB", "a_nBr"))]

corr_pearson <- cor(data_clean_sans_non_signif[,-19], method = "pearson")
corrplot(corr_pearson,method = "number", type = "lower")

data_druglike <- subset(data_clean_sans_non_signif, drug == "druglike")
data_nondruglike <- subset(data_clean_sans_non_signif, drug == "nondruglike")

druglike_indices <- sample(1:nrow(data_druglike), size = 2/3 * nrow(data_druglike))
nondruglike_indices <- sample(1:nrow(data_nondruglike), size = 2/3 * nrow(data_nondruglike))

train_druglike <- data_druglike[druglike_indices, ]
train_nondruglike <- data_nondruglike[nondruglike_indices, ]

validation_druglike <- data_druglike[-druglike_indices, ]
validation_nondruglike <- data_nondruglike[-nondruglike_indices, ]

train_dataset <- rbind(train_druglike, train_nondruglike)
validation_dataset <- rbind(validation_druglike, validation_nondruglike)

variables_quantitatives <- names(train_dataset)[sapply(train_dataset, is.numeric) & names(train_dataset) != "drug"]
vec <- c()

for(var in variables_quantitatives){
  data_app <- train_dataset[[var]]
  data_val <- validation_dataset[[var]]

  test_result <- t.test(data_app, data_val)
  vec <- c(vec, test_result$p.value)
}

names(vec) <- variables_quantitatives
print(vec)

# TOUS >0.05 donc échantillons homogènes.

par(mfrow=c(4, 5), mar=c(2, 2, 2, 2))

for(var in variables_quantitatives) {
  # Calculer la plage commune pour chaque variable pour standardiser les axes des histogrammes
  common_range <- range(c(train_dataset[[var]], validation_dataset[[var]]), na.rm = TRUE)
  
  # Histogramme pour l'échantillon d'apprentissage
  hist(train_dataset[[var]], breaks = 20, col = rgb(1, 0, 0, 0.5), xlim = common_range,
       main = paste("Histogramme de", var), xlab = var, ylab = "Fréquence")
  
  # Superposer l'histogramme pour l'échantillon de validation
  hist(validation_dataset[[var]], breaks = 20, col = rgb(0, 0, 1, 0.5), add = TRUE)
  

}

# Réinitialiser les paramètres graphiques à un seul plot par page
par(mfrow=c(1, 1))


table_app <- table(train_dataset$drug)
table_val <- table(validation_dataset$drug)
contingency_table <- rbind(Apprentissage=table_app, Validation=table_val)

chi2_result <- chisq.test(contingency_table, correct = F)
print(chi2_result)

train_dataset$drug <- factor(train_dataset$drug, levels = c("nondruglike", "druglike"))

model <- glm(drug ~ ., data = train_dataset, family = binomial)
summary(model)


# Encoder la variable cible 'drug' en 0 et 1
y_reelle <- ifelse(validation_dataset$drug == "druglike", 1, 0)

# Vérifiez si la conversion a bien fonctionné
table(y_reelle)  # Devrait montrer le nombre de 0 et de 1

# Prédire les probabilités sur l'échantillon de validation
y_pred_prob <- predict(model, newdata=validation_dataset, type="response")

# Définir le seuil pour la classification
seuil <- 0.5

# Convertir les probabilités en prédictions binaires
y_predite <- ifelse(y_pred_prob > seuil, 1, 0)


# Calculer le nombre de vrais positifs, vrais négatifs, faux positifs, et faux négatifs
nb_VP <- sum((y_predite == 1) & (y_reelle == 1))
nb_VN <- sum((y_predite == 0) & (y_reelle == 0))
nb_FP <- sum((y_predite == 1) & (y_reelle == 0))
nb_FN <- sum((y_predite == 0) & (y_reelle == 1))


taux_erreur <- (nb_FP + nb_FN) / length(y_reelle)
taux_bien_predit <- (nb_VP + nb_VN) / length(y_reelle)
sensibilite <- nb_VP / (nb_VP + nb_FN)
specificite <- nb_VN / (nb_VN + nb_FP)

# Afficher les résultats
cat("Taux d'erreur:", taux_erreur, "\n",
    "Taux bien prédit:", taux_bien_predit, "\n",
    "Sensibilité (Taux de vrais positifs):", sensibilite, "\n",
    "Spécificité (Taux de vrais négatifs):", specificite, "\n")

# Opt

# Création du modèle sans certaines variables
reduced.model <- glm(drug ~ ., data = train_dataset, family = binomial)

# Utilisation de la méthode stepwise pour la sélection de variables
step.model <- step(reduced.model, direction = "both")

# Affichage du résumé du modèle final
summary(step.model)

y_reelle <- ifelse(validation_dataset$drug == "druglike", 1, 0)

# Vérifiez si la conversion a bien fonctionné
table(y_reelle)  # Devrait montrer le nombre de 0 et de 1

# Prédire les probabilités sur l'échantillon de validation
y_pred_prob <- predict(step.model, newdata=validation_dataset, type="response")

# Définir le seuil pour la classification
seuil <- 0.5

# Convertir les probabilités en prédictions binaires
y_predite <- ifelse(y_pred_prob > seuil, 1, 0)


# Calculer le nombre de vrais positifs, vrais négatifs, faux positifs, et faux négatifs
nb_VP <- sum((y_predite == 1) & (y_reelle == 1))
nb_VN <- sum((y_predite == 0) & (y_reelle == 0))
nb_FP <- sum((y_predite == 1) & (y_reelle == 0))
nb_FN <- sum((y_predite == 0) & (y_reelle == 1))


taux_erreur <- (nb_FP + nb_FN) / length(y_reelle)
taux_bien_predit <- (nb_VP + nb_VN) / length(y_reelle)
sensibilite <- nb_VP / (nb_VP + nb_FN)
specificite <- nb_VN / (nb_VN + nb_FP)

# Afficher les résultats
cat("Taux d'erreur:", taux_erreur, "\n",
    "Taux bien prédit:", taux_bien_predit, "\n",
    "Sensibilité (Taux de vrais positifs):", sensibilite, "\n",
    "Spécificité (Taux de vrais négatifs):", specificite, "\n")


# Nom des variables dans le modèle initial
variables_initiales <- names(coef(reduced.model))

# Nom des variables dans le modèle final (après sélection)
variables_finales <- names(coef(step.model))

# Trouver les variables retirées
variables_retirees <- setdiff(variables_initiales, variables_finales)

# Afficher les variables retirées
print(variables_retirees)


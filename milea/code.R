library(car)
library(foreign)

milea = read.spss(file.choose(), to.data.frame = T)
milea1 = milea[!milea$Protein_mutation_gene_1 %in% "not applicable (NA)", c("Sex", "Protein_mutation_gene_1", "AXIN2_222696_at", "Tumor_localisation", "Tumor_size", "Recurrence", "Age")]
milea1$Protein_mutation_gene_1 = droplevels(milea1$Protein_mutation_gene_1)
milea1_model = lm(data=milea1, AXIN2_222696_at~Sex + Protein_mutation_gene_1 + Age + Tumor_localisation + Tumor_size + Recurrence)
shapiro.test(residuals(milea1_model))
plot(milea1_model)
linearHypothesis(milea1_model, matchCoefs(milea1_model, "Protein_mutation_gene_1"))


milea2_model = lm(data=milea1, AXIN2_222696_at~Protein_mutation_gene_1  + Tumor_localisation + Age + Recurrence)
plot(milea2_model)
linearHypothesis(milea2_model, matchCoefs(milea2_model, "Protein_mutation_gene_1"))


milea3_model = lm(data=milea1, AXIN2_222696_at~Protein_mutation_gene_1 + Tumor_localisation)
plot(milea3_model)
linearHypothesis(milea3_model, matchCoefs(milea3_model, "Protein_mutation_gene_1"))

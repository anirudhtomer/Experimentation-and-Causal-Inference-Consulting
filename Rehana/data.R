library(foreign)
#Read the whole dataset
rahman = read.spss(file=file.choose(), to.data.frame = T)
rahman_under = droplevels(rahman[!is.na(rahman$Reduction_UN),])
rahman_over = droplevels(rahman[!is.na(rahman$Reduction_OV),])

by(data=rahman$Reduction_UN_OR_OV, INDICES = rahman$Group, function(x){table(x)})

summary(by(data=rahman[rahman$Group=="controle",], INDICES = rahman[rahman$Group=="controle",]$Subjectnr, nrow))
summary(by(data=rahman[rahman$Group!="controle",], INDICES = rahman[rahman$Group!="controle",]$Subjectnr, nrow))


by(data=rahman$Reduction_UN_OR_OV, INDICES = rahman$Group, function(x){table(x)/length(x)})

glmer(Reduction_UN_OR_OV~Group + (1|Subjectnr), data=rahman, family=binomial,nAGQ=20)

dim_mod = mixed_model(fixed = Reduction_UN_OR_OV ~ Group, 
            random = ~ 1 | Subjectnr, data = rahman, 
            family = binomial())

marginal_dim_mod = marginal_coefs(dim_mod, std_errors = TRUE)
marginal_dim_mod

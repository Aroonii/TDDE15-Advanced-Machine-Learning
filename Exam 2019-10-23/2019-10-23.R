library(kernlab)
library(gRain)
library(bnlearn)

monty_network = model2network("[D][P][M|D:P]")



custom.fit(monty_network, data = "A", dist = names(monty_network))
cptD = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("D1", "D2", "D3")))
cptP = matrix(c(1/3, 1/3, 1/3), ncol = 3, dimnames = list(NULL, c("P1", "P2", "P3")))
cptM = c(
  0,.5,.5,
  0,0,1,
  0,1,0,
  0,0,1,
  .5,0,.5,
  1,0,0,
  0,1,0,
  1,0,0,
  .5,.5,0)
dim(cptM) = c(3,3,3)
dimnames(cptM) = list("M" = c("M1", "M2", "M3"), "D" =  c("D1", "D2", "D3"), "P" = c("P1", "P2", "P3"))
MHfit = custom.fit(monty_network, list(D = cptD, P = cptP, M = cptM))
MHcom = compile(as.grain(MHfit))
plot(MHcom)

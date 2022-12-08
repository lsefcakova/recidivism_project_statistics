library(ggplot2)
library(latex2exp)



dfdbetaj <- function(x,lambda,beta_MLE){
  return(x + lambda*x/abs(x) - beta_MLE)
}

ggplot() + 
  xlim(-1, 5) + 
  geom_function(fun=dfdbetaj, args = list(lambda = 1, beta_MLE = 2), col = 'blue') +
  geom_function(fun=dfdbetaj, args = list(lambda = 3, beta_MLE = 2), col = 'red') +
  xlab('beta_j') +
  geom_hline(aes(yintercept = 0), linetype =2) +
  geom_vline(aes(xintercept = 2) , linetype =2) +
  xlab(TeX("$\\beta_j$")) +
  ylab(TeX("$\\frac{\\partial f}{\\partial \\beta_j}$")) +
  annotate("text", x = 4, y = 5.9, label = TeX("$\\lambda = 3$"), size = 3) +
  annotate("text", x = 4, y = 2.3, label = TeX("$\\lambda = 1$"), size = 3) +
  annotate("text", x = 2.3, y = -5, label = TeX("$\\beta _{MLE}$"), size = 3) +
  ggtitle(TeX("Case $\\beta _{MLE} > 0 $")) +
  theme(axis.title.y = element_text(angle=0)) +
  ggsave('derivative_lasso_funct_OLSpos.png')


ggplot() + 
  xlim(-5, 1) + 
  geom_function(fun=dfdbetaj, args = list(lambda = 1, beta_MLE = -2), col = 'blue') +
  geom_function(fun=dfdbetaj, args = list(lambda = 3, beta_MLE = -2), col = 'red') +
  xlab('beta_j') +
  geom_hline(aes(yintercept = 0), linetype =2) +
  geom_vline(aes(xintercept = -2) , linetype =2) +
  xlab(TeX("$\\beta_j$")) +
  ylab(TeX("$\\frac{\\partial f}{\\partial \\beta_j}$")) +
  annotate("text", x = -3.5, y = -1.5, label = TeX("$\\lambda = 1$"), size = 3) +
  annotate("text", x = -3.5, y = -5.5, label = TeX("$\\lambda = 3$"), size = 3) +
  annotate("text", x = -2.3, y = -5, label = TeX("$\\beta _{MLE}$"), size = 3) +
  ggtitle(TeX("Case $\\beta _{MLE} < 0 $")) +
  theme(axis.title.y = element_text(angle=0)) +
  ggsave('derivative_lasso_funct_OLSneg.png')


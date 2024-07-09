# s2 Intuition: ----------------------------------------------------------------
nTaxa <- 10
n_iter <- 1000
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)

s2Intui <- lapply(1:length(s2List), function(y){
  sapply(1:n_iter, function(x){
    set.seed(x + 1)
    xiDum <- rnorm(nTaxa, 0, sqrt(s2List[y]))
    relativeDum <- exp(xiDum)/sum(exp(xiDum))
    var(exp(xiDum)/sum(exp(xiDum)))
  }) %>% as.data.frame()
}) %>%
  bind_cols() %>%
  mutate(Iter = 1:n_iter) %>%
  pivot_longer(!Iter)

s2Intui$name <- factor(s2Intui$name, levels = paste0("....", 1:length(s2List)),
                       labels = paste0("s2 = ", s2List))

s2Intui %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot()

nTaxaList <- c(10, 25, 50, 100, 250, 500)
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)
nS2 <- expand.grid(n = nTaxaList, s2 = s2List)
n_iter <- 1000

s2Intui <- lapply(1:nrow(nS2), function(y){
  
  lapply(1:n_iter, function(x){
    set.seed(x)
    xiDum <- rnorm(nS2[y, 1], 0, sqrt(nS2[y, 2]))
    relativeDum <- exp(xiDum)/sum(exp(xiDum))
    c(Iter = x, n = nS2[y, 1], s2 = nS2[y, 2], var = var(relativeDum))
  }) %>%
    bind_rows()
  
}) %>%
  bind_rows() %>%
  transmute(n = paste0("n = ", n), s2 = paste0("s2 = ", s2), var)

s2Intui$n <- factor(s2Intui$n, levels = paste0("n = ", nTaxaList))
s2Intui$s2 <- factor(s2Intui$s2, levels = paste0("s2 = ", s2List))

ggplot(s2Intui, aes(x = s2, y = var)) +
  geom_boxplot() +
  facet_wrap(. ~ n, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 270))

nTaxa <- 50
s2List <- c(1e-3, 0.01, 0.1, 1, 2.5, 5, 10, 25, 50, 100)
set.seed(1)
x <- rnorm(nTaxa, 0, sqrt(1))

data.frame(j = 1:50, p = exp(x)/sum(exp(x))) %>%
  ggplot(aes(x = factor(j), y = p)) +
  geom_bar(stat = "identity")


### n_xi intuition: ------------------------------------------------------------
data.frame(x = paste0("OTU ", 1:60),
           actual_y = as.numeric(colSums(otuHIV)/sum(otuHIV))) %>%
  ggplot(aes(x = x, y = actual_y)) +
  geom_bar(stat = "identity")

x <- rnorm(60, 0, sqrt(1))
estimated_Y <- exp(x)/sum(exp(x))

var(as.numeric(colSums(otuHIV)/sum(otuHIV)))

data.frame(x = paste0("OTU ", 1:60), estimated_Y) %>%
  ggplot(aes(x = x, y = estimated_Y)) +
  geom_bar(stat = "identity")


n <- 60
s2 <- 1
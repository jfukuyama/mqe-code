library(ape)
library(magrittr)
library(phyloseq)
library(latex2exp)
source("eigen-fns.R")
########## balanced binary tree ##########
bbt = stree(n = 32, type = "balanced")
bbt$edge.length = rep(1, nrow(bbt$edge))

## + 1 at the end because we want some length associated with the root
Qbbt = vcv(bbt) + 1
Qbbt_eig = eigen(Qbbt)

## the eigenvectors you get this way are linear combinations of the more interpretable eigenvectors
df = data.frame(evec = as.vector(Qbbt_eig$vector),
                tip_index = rep(1:nrow(Qbbt_eig$vector), times = ncol(Qbbt_eig$vector)),
                evec_index = rep(1:ncol(Qbbt_eig$vector), each = nrow(Qbbt_eig$vector)))

p1 = ggplot(subset(df, evec_index  %in% c(1,2,3,4,5, 31, 32)))+
    geom_point(aes(x = tip_index, y = evec)) +
    facet_grid(evec_index ~ .) +
    scale_y_continuous(breaks = c(-.5, 0, .5), limits = c(-.8, .8), name = "")+
    scale_x_continuous(name="")
tree_plot = phyloseq::plot_tree(bbt) + coord_flip() + scale_x_reverse()

## so we plot the theoretical ones instead
type_df = data.frame(i = c(5,4,3,3,1,1,1,0), j = 5, k = c(0,0,0,1,0,1,2,15))
type_orders = mapply(get_evec_order, type_df$i, type_df$j, type_df$k)
evecs = mapply(cijk, type_df$i, type_df$j, type_df$k) %>% as.vector
df = data.frame(
    x = rep(1:32, times = nrow(type_df)),
    y = evecs,
    type = rep(type_orders, each = 32))


p1 = ggplot(df)+ geom_point(aes(x = x, y = y)) +
    facet_grid(type ~ .) +
    scale_y_continuous(breaks = c(-1, 0, 1), limits = c(-1.2, 1.2), name = "") +
    scale_x_continuous(name = "")
p_whole = treeDA::combine_plot_and_tree(p1, tree_plot, tree.height = 3)
plot(p_whole)


r_vals = c(0,.1, .5, .95,  1)
evals = lapply(r_vals, get_evals, Qbbt_eig$values)
eval_df_bbt = mapply(FUN = function(evals, r) data.frame(evals = evals / sum(evals), idx = 1:length(evals), r = paste0("r = ", r)), evals, r_vals, SIMPLIFY = FALSE) %>% Reduce(rbind, .)
p_evals = ggplot(eval_df_bbt) +
    geom_point(aes(x = idx, y = evals)) +
    facet_grid(. ~ r) +
    scale_x_continuous("") + scale_y_continuous("Fraction of trace")

plot(p_evals)

########## comb tree ##########
p = 32
ct = stree(n = p, type = "left")
ct$edge.length = rep(1, nrow(ct$edge)) / p

Qct = vcv(ct)
Qct_eig = eigen(Qct)

df = data.frame(evec = as.vector(Qct_eig$vector),
                tip_index = rep(1:nrow(Qct_eig$vector), times = ncol(Qct_eig$vector)),
                evec_index = rep(1:ncol(Qct_eig$vector), each = nrow(Qct_eig$vector)))

p1 = ggplot(subset(df, evec_index  %in% c(1:8)))+
    geom_point(aes(x = tip_index, y = evec)) +
    facet_grid(evec_index ~ .) +
    scale_y_continuous(breaks = c(-.5, 0, .5), name = "")+
    scale_x_continuous(name="")
tree_plot = phyloseq::plot_tree(ct) + coord_flip() + scale_x_reverse()

p_whole = treeDA::combine_plot_and_tree(p1, tree_plot, tree.height=3)
plot(p_whole)
r_vals = c(0,.1, .5, .95,  1)
evals = lapply(r_vals, get_evals, Qct_eig$values)
eval_df_ct = mapply(FUN = function(evals, r) data.frame(evals = evals / sum(evals), idx = 1:length(evals), r = paste0("r = ", r)), evals, r_vals, SIMPLIFY = FALSE) %>% Reduce(rbind, .)
p_evals = ggplot(eval_df_ct) +
    geom_point(aes(x = idx, y = evals)) +
    facet_grid(. ~ r) +
    scale_x_continuous("") + scale_y_continuous("Fraction of trace")

plot(p_evals)


########## real tree ##########

data(AntibioticPhyloseq, package = "adaptiveGPCA")
Qrt = vcv(phy_tree(AntibioticPhyloseq))
Qrt_eig = eigen(Qrt)
max_vals = apply(Qrt_eig$vectors, 2, function(x) max(abs(x)))
Qrt_evecs_scaled = sweep(Qrt_eig$vectors, MARGIN = 2, STATS = max_vals, FUN = '/')
df = data.frame(evec = as.vector(Qrt_evecs_scaled),
                tip_index = rep(1:nrow(Qrt_eig$vector), times = ncol(Qrt_eig$vector)),
                evec_index = rep(1:ncol(Qrt_eig$vector), each = nrow(Qrt_eig$vector)))

p1 = ggplot(subset(df, evec_index  %in% c(1,2,4,5,50,70,500,1500)))+
    geom_point(aes(x = tip_index, y = evec)) +
    facet_grid(evec_index ~ .) +
    scale_y_continuous(breaks = c(-.5, 0, .5), name = "")+
    scale_x_continuous(name="")
tree_plot = phyloseq::plot_tree(phy_tree(AntibioticPhyloseq)) + coord_flip() + scale_x_reverse()
p_whole = treeDA::combine_plot_and_tree(p1, tree_plot, tree.height=3)

########## eigenvalue plot ##########

eigen_list = list("Balanced binary" = Qbbt_eig, "Comb tree" = Qct_eig, "Real tree" = Qrt_eig)
eval_info_list = lapply(eigen_list, function(x) {
    data.frame(eval = x$values / sum(x$values), index = 1:length(x$values))
})
eval_info_list = mapply(FUN = function(x, y) data.frame(x, type = y),
                        eval_info_list, names(eigen_list), SIMPLIFY = FALSE)
eval_df = Reduce(rbind, eval_info_list)
ggplot(eval_df) + geom_point(aes(x = index, y = eval)) +
    facet_wrap( ~ type, scales = "free_x", ncol = 2) +
    scale_x_continuous("") + scale_y_continuous("Normalized eigenvalue")

########## big outgroup ##########
bot = stree(n = 100, type = "left")
bot$edge.length = c(rep(1, nrow(bot$edge) - 1), 300)

Qbot = vcv(bot)
Qbot_eig = eigen(Qbot)

df = data.frame(evec = as.vector(Qbot_eig$vector),
                tip_index = rep(1:nrow(Qbot_eig$vector), times = ncol(Qbot_eig$vector)),
                evec_index = rep(1:ncol(Qbot_eig$vector), each = nrow(Qbot_eig$vector)))

p1 = ggplot(subset(df, evec_index  %in% c(1,2,3,4,5, 50,98, 99, 100)))+
    geom_point(aes(x = tip_index, y = evec)) +
    facet_grid(evec_index ~ .) +
    scale_y_continuous(breaks = c(-.5, 0, .5), name = "")+
    scale_x_continuous(name="")
tree_plot = phyloseq::plot_tree(bot) + coord_flip() + scale_x_reverse()

p_whole = treeDA::combine_plot_and_tree(p1, tree_plot, tree.height=3)
plot(p_whole)

r_vals = c(0,.1, .5, .95,  1)
evals = lapply(r_vals, get_evals, Qbot_eig$values)
eval_df_bot = mapply(FUN = funbotion(evals, r) data.frame(evals = evals / sum(evals), idx = 1:length(evals), r = paste0("r = ", r)), evals, r_vals, SIMPLIFY = FALSE) %>% Reduce(rbind, .)
p_evals = ggplot(eval_df_bot) +
    geom_point(aes(x = idx, y = evals)) +
    facet_grid(. ~ r) +
    scale_x_continuous("") + scale_y_continuous("Frabotion of trace")

plot(p_evals)

########## Checking on the eigenvalues of the balanced binary tree ##########
bbt = stree(n = 32, type = "balanced")
bbt$edge.length = rep(1, nrow(bbt$edge))

Qbbt = vcv(bbt)
Qbbt_eig = eigen(Qbbt)




T = Dij(0,5) + Dij(1,5) + Dij(2,5) + Dij(3, 5) + Dij(4,5)
sapply(0:4, function(i) sum(2^(0:i))) %>% sort
unique(eigen(T)$values) %>% sort

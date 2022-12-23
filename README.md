# mdtm: Multiple Descriptors Topic Models

## Installation
To install this package type (in R):

```
devtools::install_github("giovannitoto/GGIF")
```

## Contents
The package allows to run a Collapsed Gibbs Sampler (CGS) to perform the posterior inference of the following topic models:
* [Latent Dirichlet Allocation (LDA)](https://dl.acm.org/doi/10.5555/944919.944937)
* [Twitter-LDA (TLDA)](https://link.springer.com/chapter/10.1007/978-3-642-20161-5_34)
* [Hashtag-LDA (HLDA)](https://www.sciencedirect.com/science/article/abs/pii/S0167739X15003258?via%3Dihub)
* [Microblog-LDA, MLDA](https://ceur-ws.org/Vol-3177/paper1.pdf)

## Example
We want to estimate *LDA*, *TLDA*, *HLDA* and *MLDA* on a collection of documents and then predict the latent structure of new documents.

We set the number of topics to $T=30$ and we consider symmetric Dirichlet prior distributions with parameter $50/T$ for the user-topic distributions $\theta^*$ and the document-topic distributions $\theta$. Then, since this is a simple example, we set the number of iterations to 100 and the seed to 28.
```
TOPICS <- 30
alpha <- alphastar <- rep(1, TOPICS) * 50 / TOPICS
iterations <- 100
seed <- 28
```

### Inference via CGS
We launch the CGSs to perform the posterior inference of the four topic models: the first arguments are related to the models, i.e. are the parameters of the prior distributions, while the remaining ones are related to algorithm. All the following functions work in the same way that is:
1. create a folder, called `result_folder`, in which to save all the files mentioned below;
2. create a `hyperparameter.RDS` file containing the arguments of the function and other useful quantities;
3. for each iteration create a subfolder containing the state of the chain at that iteration.

```
mdtm::CGS_LDA(w = train$t_matrix, alpha = alpha, betaV = rep(0.1, max(train$t_matrix)),
              iterations = iterations, seed = seed,
              result_folder = file.path("results", "libLDA"))
mdtm::CGS_TwitterLDA(w = train$t_matrix, doc_users = train$doc_users,
                     alphastar = alphastar, betaV = rep(0.1, max(train$t_matrix)),
                     bV = c(1, 1), iterations = iterations, seed = seed,
                     result_folder = file.path("results", "libTLDA"))
mdtm::CGS_HashtagLDA(w = train$w_matrix, h = train$h_matrix, doc_users = train$doc_users,
                     alphastar = alphastar, betaV = rep(0.1, max(train$w_matrix)),
                     betaH = rep(0.1, max(train$h_matrix)), bH = c(1, 1),
                     iterations = iterations, seed = seed,
                     result_folder = file.path("results", "libHLDA"))
mdtm::CGS_MicroblogLDA(w = list(train$w_matrix, train$h_matrix), doc_users = train$doc_users,
                       alphastar = alphastar, alpha = alpha,
                       beta = list(rep(0.1, max(train$w_matrix)), rep(0.1, max(train$h_matrix))),
                       b = list(c(1, 1), c(1, 1)), bdelta = c(1, 1), bT = c(1, 1),
                       alpha0 = 10^-7, iterations = iterations, seed = seed,
                       result_folder = file.path("results", "libMLDA"))
```

### Posterior estimates
Now we have the chains of the latent variables, however we do not have the posterior estimates of the parameters, such as $\theta^*,\theta,\ldots$. As before, all the following functions work in the same way that is:
1. compute the logarithm of the at each iteration;
2. compute the posterior estimates of the parameters and of the latent variables;
3. save the results in a RDS file called `postproc_file.RDS`.

```
mdtm::postproc_LDA(result_folder = file.path("results", "libLDA"),
                   postproc_file = "postproc", iterations = NULL, verbose = TRUE)

mdtm::postproc_TwitterLDA(result_folder = file.path("results", "libTLDA"),
                          postproc_file = "postproc", iterations = NULL, verbose = TRUE)

mdtm::postproc_HashtagLDA(result_folder = file.path("results", "libHLDA"),
                          postproc_file = "postproc", iterations = NULL, verbose = TRUE)

mdtm::postproc_MicroblogLDA(result_folder = file.path("results", "libMLDA"),
                            postproc_file = "postproc", iterations = NULL, verbose = TRUE)
```

### Prediction via CGS
Now that we have the chains and the posterior estimates of the latent structure, we can predict the latent structure of new documents. Two approaches are available:
- Independently resample topics for new documents
```
mdtm::pred_LDA(w = test$t_matrix, betaV_new = rep(0.1, max(test$t_matrix)-max(train$t_matrix)),
               postproc_file = file.path("results", "libLDA", "postproc.RDS"),
               single_doc = TRUE, iterations = 20, seed = 28,
               result_folder = file.path("results", "predLDA_single"))
mdtm::pred_TwitterLDA(w = test$t_matrix, doc_users = test$doc_users,
                      betaV_new = rep(0.1, max(test$t_matrix)-max(train$t_matrix)),
                      postproc_file = file.path("results", "libTLDA", "postproc.RDS"),
                      single_doc = TRUE, iterations = 20, seed = 28,
                      result_folder = file.path("results", "predTLDA_single"))
mdtm::pred_HashtagLDA(w = test$w_matrix, h = test$h_matrix, test$doc_users,
                      betaV_new = rep(0.1, max(test$w_matrix)-max(train$w_matrix)),
                      betaH_new = rep(0.1, max(test$h_matrix)-max(train$h_matrix)),
                      postproc_file = file.path("results", "libHLDA", "postproc.RDS"),
                      single_doc = TRUE, iterations = 20, seed = 28,
                      result_folder = file.path("results", "predHLDA_single"))
mdtm::pred_MicroblogLDA(w = list(test$w_matrix, test$h_matrix), doc_users = test$doc_users,
                        beta_new = list(rep(0.1, max(test$w_matrix)-max(train$w_matrix)),
                                        rep(0.1, max(test$h_matrix)-max(train$h_matrix))),
                        postproc_file = file.path("results", "libMLDA", "postproc.RDS"),
                        single_doc = TRUE, iterations = 20, seed = 28,
                        result_folder = file.path("results", "predMLDA_single"))
```
- Jointly resample topics for all new documents
```
mdtm::pred_LDA(w = test$t_matrix, betaV_new = rep(0.1, max(test$t_matrix)-max(train$t_matrix)),
               postproc_file = file.path("results", "libLDA", "postproc.RDS"),
               single_doc = FALSE, iterations = 20, seed = 28,
               result_folder = file.path("results", "predLDA_all"))
mdtm::pred_TwitterLDA(w = test$t_matrix, doc_users = test$doc_users,
                      betaV_new = rep(0.1, max(test$t_matrix)-max(train$t_matrix)),
                      postproc_file = file.path("results", "libTLDA", "postproc.RDS"),
                      single_doc = FALSE, iterations = 20, seed = 28,
                      result_folder = file.path("results", "predTLDA_all"))
mdtm::pred_HashtagLDA(w = test$w_matrix, h = test$h_matrix, test$doc_users,
                      betaV_new = rep(0.1, max(test$w_matrix)-max(train$w_matrix)),
                      betaH_new = rep(0.1, max(test$h_matrix)-max(train$h_matrix)),
                      postproc_file = file.path("results", "libHLDA", "postproc.RDS"),
                      single_doc = FALSE, iterations = 20, seed = 28,
                      result_folder = file.path("results", "predHLDA_all"))
mdtm::pred_MicroblogLDA(w = list(test$w_matrix, test$h_matrix), doc_users = test$doc_users,
                        beta_new = list(rep(0.1, max(test$w_matrix)-max(train$w_matrix)),
                                        rep(0.1, max(test$h_matrix)-max(train$h_matrix))),
                        postproc_file = file.path("results", "libMLDA", "postproc.RDS"),
                        single_doc = FALSE, iterations = 20, seed = 28,
                        result_folder = file.path("results", "predMLDA_all"))
```
In the first approach, the latent structure of a new document depends on the latent structure of the known documents; in the second approach, the latent structure of a new document depends on the latent structure of both known and new documents.


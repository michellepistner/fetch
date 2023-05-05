# fetch
Welcome to the `fetch` repository on GitHub!

`fetch` is built on the `ALDEx2` package but also allows for global test statistics to be computed. These test statistics can be used to determine which features should be tested on a taxa by taxa basis.

To run `fetch` with one condition:

```
library(fetch)
data(selex)
#subset for efficiency
selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))
f <- fetch(selex, conds, mc.samples=10, denom="all", test = "global")
f
```

To run `fetch` with more than one condition:

```
data(selex)
#subset for efficiency
selex <- selex[1201:1600,]
covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
                         "B" = c(rep(0, 7), rep(1, 7)))
mm <- model.matrix(~ A + B, covariates)
f <- fetch(selex, mm, mc.samples=10, denom="all", test = "global")
f
```

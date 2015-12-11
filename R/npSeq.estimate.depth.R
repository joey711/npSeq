npSeq.estimate.depth <- function(n)
{
	iter <- 5
	cmeans <- colSums(n) / sum(n)
	for (i in 1:iter)
	{
		n0 <- rowSums(n) %*% t(cmeans)
		prop <- rowSums((n - n0)^2 / (n0 + 1e-08))
		qs <- quantile(prop, c(0.25, 0.75))
		keep <- (prop >= qs[1]) & (prop <= qs[2])
		cmeans <- colMeans(n[keep, ])
		cmeans <- cmeans / sum(cmeans)
	}
	depth <- cmeans / mean(cmeans)
	
	return(depth)
}

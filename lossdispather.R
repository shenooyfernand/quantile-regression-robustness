loss_dispatch <- function(loss, r, X, tau, k = NULL, c = NULL) {
  switch(loss,
         "pinball" = list(
           loss = pinball_loss(r, tau),
           grad = pinball_grad(r, X, tau)
         ),
         "qhuber" = list(
           loss = qhuber_loss(r, tau, k),
           grad = qhuber_grad(r, X, tau, k)
         ),
         "logcosh" = list(
           loss = logcosh_loss(r, tau, c),
           grad = logcosh_grad(r, X, tau, c)
         )
  )
}

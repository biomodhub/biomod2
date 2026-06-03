summary(ur.df(Canada[, "prod"],
              type = "trend", lags = 2))
summary(ur.df(diff(Canada[, "prod"]),
              type = "drift", lags = 1))
summary(ur.df(Canada[, "e"],
              type = "trend", lags = 2))
summary(ur.df(diff(Canada[, "e"]),
              type = "drift", lags = 1))
summary(ur.df(Canada[, "U"],
              type = "drift", lags = 1))
summary(ur.df(diff(Canada[, "U"]),
              type = "none", lags = 0))
summary(ur.df(Canada[, "rw"],
              type = "trend", lags = 4))
summary(ur.df(diff(Canada[, "rw"]),
              type = "drift", lags = 3))
summary(ur.df(diff(Canada[, "rw"]),
              type = "drift", lags = 0))


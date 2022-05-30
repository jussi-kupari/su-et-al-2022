nuke_from_orbit <- function() {
  message("It's the only way to be sure.")

  kaboom <- r"{


                                  ________________
                             ____/ (  (    )   )  \___
                            /( (  (  )   _    ))  )   )\
                          ((     (   )(    )  )   (   )  )
                        ((/  ( _(   )   (   _) ) (  () )  )
                       ( (  ( (_)   ((    (   )  .((_ ) .  )_
                      ( (  )    (      (  )    )   ) . ) (   )
                     (  (   (  (   ) (  _  ( _) ).  ) . ) ) ( )
                     ( (  (   ) (  )   (  ))     ) _)(   )  )  )
                    ( (  ( \ ) (    (_  ( ) ( )  )   ) )  )) ( )
                     (  (   (  (   (_ ( ) ( _    )  ) (  )  )   )
                    ( (  ( (  (  )     (_  )  ) )  _)   ) _( ( )
                     ((  (   )(    (     _    )   _) _(_ (  (_ )
                      (_((__(_(__(( ( ( |  ) ) ) )_))__))_)___)
                      ((__)        \\||lll|l||///          \_))
                               (   /(/ (  )  ) )\   )
                             (    ( ( ( | | ) ) )\   )
                              (   /(| / ( )) ) ) )) )
                            (     ( ((((_(|)_)))))     )
                             (      ||\(|(|)|/||     )
                           (        |(||(||)||||        )
                             (     //|/l|||)|\\ \     )
                           (/ / //  /|//||||\\  \ \  \ _)

    888    d8P         d8888 888888b.    .d88888b.   .d88888b.  888b     d888 888
    888   d8P         d88888 888  "88b  d88P" "Y88b d88P" "Y88b 8888b   d8888 888
    888  d8P         d88P888 888  .88P  888     888 888     888 88888b.d88888 888
    888d88K         d88P 888 8888888K.  888     888 888     888 888Y88888P888 888
    8888888b       d88P  888 888  "Y88b 888     888 888     888 888 Y888P 888 888
    888  Y88b     d88P   888 888    888 888     888 888     888 888  Y8P  888 Y8P
    888   Y88b   d8888888888 888   d88P Y88b. .d88P Y88b. .d88P 888   "   888  "
    888    Y88b d88P     888 8888888P"   "Y88888P"   "Y88888P"  888       888 888




  }"

  cat(kaboom)

  pacs <- names(sessionInfo()$otherPkgs)

  if (!(is.null(pacs))) {
    suppressWarnings(invisible(lapply(
      paste("package:", pacs, sep = ""),
      detach,
      character.only = TRUE,
      unload = TRUE
    )))
  }

  suppressWarnings(rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv))

  .rs.restartR()
}

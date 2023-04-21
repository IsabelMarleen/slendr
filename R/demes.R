#' Load and validate a Demes model in slendr format
#'
#' @param demes_path File path to yaml file specifying the Demes model
#'
#' @return Compiled \code{slendr_model} model object which encapsulates all
#'   information about the specified model (which populations are involved,
#'   when and how much gene flow should occur, what is the spatial resolution
#'   of a map, and what spatial dispersal and mating parameters should be used
#'   in a SLiM simulation, if applicable)
#' @export
#'
#' @examples
#' path_jacobs <- system.file("extdata", "models", "Demes", "jacobs.yaml", package = "slendr")
#' slendr_model <- compile_demes(path_jacobs)
compile_demes <- function(demes_path, output_path = NULL){
  demes <- demes::read_demes(demes_path)

  # TODO: Maybe introduce checks for things like infinite start_times where demes and slendr disagree

  # Populations
  # TODO: Fix Parents
  pops <- list()
  anchor <- find_oldest_time_anchor(demes)+1
  gen_time <- demes$generation_time
  for (i in 1:length(demes$demes)){
    pops[[i]] <- convert_deme(demes$demes[[i]], pops, anchor)
    names(pops)[i] <- pops[[i]]$pop
  }

  # Resize / Epochs
  pops <- purrr::map(demes$demes, ~convert_epochs(.x$epochs, pop=pops[[.x$name]], gen_time))

  # Geneflow events
  gf_migs <- purrr::map(demes$migrations, convert_migration, pops=pops, gen_time)
  gf_pulses <- purrr::map(demes$pulses, convert_pulse, pops=pops, gen_time)
  gf <- c(gf_migs, gf_pulses)

  if (length(gf) == 0){
    model <- compile_model(populations = pops,
                           generation_time = demes$generation_time,
                           description = demes$description,
                           path=output_path,
                           direction = "backward")
  } else {
    model <- compile_model(populations = pops,
                           generation_time = demes$generation_time,
                           gene_flow = gf,
                           description = demes$description,
                           path=output_path,
                           direction = "backward")
  }

  return(model)
}

convert_deme <- function(deme, pops, oldest_anchor){
  p_name <- deme$name
  if (deme$start_time == Inf){
    p_time <- oldest_anchor # TODO: Extract oldest time in the model as an anchor point
  } else {
    p_time <- deme$start_time
  }
  latest_epoch <- length(deme$epochs)
  if (deme$epochs[[latest_epoch]]$end_time == 0){
    p_remove <- -1
  } else {
    p_remove <- deme$epochs[[latest_epoch]]$end_time
  }

  p_start_size <- deme$epochs[[1]]$start_size

  if (length(deme$ancestors) != 0){
    p_parent <- pops[[deme$ancestors[1]]]
    pop <- population(name = p_name, time = p_time, N = p_start_size, parent = p_parent, remove=p_remove)
   }
  else {
    pop <- population(name = p_name, time = p_time, N = p_start_size, remove=p_remove)
  }

  return(pop)
}

convert_epochs <- function(epochs, pop, gen_time){
  for (i in 1:length(epochs)){
    epoch <- epochs[[i]]
    r_N <- epoch$end_size
    if (epoch$size_function == "exponential"){
      r_how <- "exponential"
    } else if (epoch$size_function == "constant" & i == 1 ){
      next
    } else {
      r_how <- "step"
    }
    r_time <- epoch$start_time
    if (r_time == Inf){
      r_time <- pop$time - 1*gen_time
    } else if (r_time == pop$time){
      r_time <- r_time - 1*gen_time
    }
    r_end <- epoch$end_time

    pop <- resize(pop, N = r_N, how = r_how, time = r_time, end = r_end)
  }
  return(pop)
}

convert_migration <- function(migration, pops, gen_time){
  if (length(migration) == 0){
    return()
  } else {
    gf_from <- pops[[migration$source]]
    gf_to <- pops[[migration$dest]]
    gf_rate <- migration$rate
    gf_start <- migration$start_time
    if (gf_start == pops[[migration$source]]$time){
      gf_start <- gf_start-(1*gen_time)
    } else if (gf_start == pops[[migration$dest]]$time){
      gf_start <- gf_start-(1*gen_time)
    }
    gf_end <- migration$end_time

    g <- gene_flow(from = gf_from, to = gf_to, rate = gf_rate, start = gf_start, end = gf_end)

    return(g)
  }
}

convert_pulse <- function(pulse, pops, gen_time){
  if (length(pulse) == 0){
    return()
  } else {
    gf_from <- pops[[pulse$source]]
    gf_to <- pops[[pulse$dest]]
    gf_rate <- pulse$proportions
    gf_start <- pulse$time
    if (gf_start == pops[[pulse$source]]$time){
      gf_start <- gf_start-(1*gen_time)
    } else if (gf_start == pops[[pulse$dest]]$time){
      gf_start <- gf_start-(1*gen_time)
    }
    gf_end <- gf_start-(1*gen_time)

    g <- gene_flow(from = gf_from, to = gf_to, rate = gf_rate, start = gf_start, end = gf_end)

    return(g)
  }
}

find_youngest_time_anchor <- function(demes){
  anchor <- min()

  return(anchor)
}

find_oldest_time_anchor <- function(demes){
  anchor_start_time <- purrr::map(demes$demes, ~.x$start_time) %>%
    unlist() %>%
    .[. != Inf]

  anchor_end_time <-  purrr::map(demes$demes, ~.x$epochs[[1]]$end_time) %>%
    unlist()

  anchor_migs <- purrr::map(demes$migrations, ~.x$start_time) %>%
    unlist()

  anchor_pulse <- purrr::map(demes$pulses, ~.x$time) %>%
    unlist()

  anchor <-  max(anchor_start_time, anchor_end_time, anchor_migs, anchor_pulse, na.rm = TRUE)

  return(anchor)
}
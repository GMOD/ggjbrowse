# ggjbrowse in grammar-of-graphics terms. Every table the engine returns is
# data, every mark kind is a geom that draws one table, and each geom paints
# the display's colors as they are, I(color), so they train no scale and the
# caller's own mappings take every scale. The keys are guides of their own
# (guides_jbrowse()). `autoplot()` only composes these pieces.

# A geom's `data` as it hands it on: a data frame is taken as the plot's
# result's `table`, so the geom draws it as it draws that table, which a geom
# drawing from several tables (`table` NULL) cannot do
replacing <- function(data, table, geom) {
  if (!is.data.frame(data)) return(data)
  if (is.null(table)) {
    stop(geom, "() draws from several of a result's tables, so its data is ",
      "a jbrowse result: change the tables in the result and pass that",
      call. = FALSE)
  }
  if (table %in% unlist(lapply(FAMILIES, names))) {
    rlang::abort(sprintf(paste0("`%s` names the records other tables key ",
      "their rows by, so a filtered copy drawn alone leaves their marks ",
      "behind: assign it into the result, `x$%s <- subset(x$%s, ...)`, ",
      "which drops its children's rows too, and draw that result"),
      table, table, table), class = "ggjbrowse_parent_data")
  }
  structure(list(frame = data, table = table), class = "jbrowse_table")
}

# What a geom draws: the jbrowse result its plot or its `data` holds, on its
# axis and read through `pick`, or that result with one table replaced by the
# data frame the caller handed the geom
table_from <- function(data, pick, finish = function(d, model) d) {
  if (inherits(data, "jbrowse_result")) {
    return(finish(pick(on_axis(data)), data))
  }
  if (is.null(data) || inherits(data, "jbrowse_table")) {
    return(function(plot_data) {
      model <- attr(plot_data, "jbrowse")
      if (is.null(model)) {
        stop("a jbrowse geom needs ggplot(<jbrowse result>) or data = ",
          "<jbrowse result>; a data frame as data replaces one of the ",
          "plot's result's tables", call. = FALSE)
      }
      if (!is.null(data)) model[[data$table]] <- data$frame
      finish(pick(on_axis(model)), model)
    })
  }
  data
}

with_defaults <- function(mapping, defaults) {
  for (name in names(mapping)) defaults[[name]] <- mapping[[name]]
  defaults
}

# A layer of a jbrowse table; one of pileup `rows` draws each at the
# display's read height
jbrowse_layer <- function(geom, data, pick, mapping, defaults, params,
                          rows = FALSE) {
  finish <- if (rows) with_inset else function(d, model) d
  layer <- ggplot2::layer(geom = geom, stat = "identity",
    position = "identity", data = table_from(data, pick, finish),
    mapping = with_defaults(mapping, defaults), inherit.aes = FALSE,
    params = params)
  ggplot2::ggproto(NULL, layer, setup_layer = function(self, data, plot) {
    warn_plot_mapping(plot$mapping)
    ggplot2::ggproto_parent(layer, self)$setup_layer(data, plot)
  })
}

# Each jbrowse geom draws its own table, so it takes no plot-level mapping;
# say so rather than let `ggplot(a, aes(fill = mapq))` draw nothing new
warn_plot_mapping <- function(mapping) {
  if (!length(mapping)) return(invisible())
  rlang::warn(paste0("jbrowse geoms draw their own tables and take no ",
    "ggplot(aes(", paste(names(mapping), collapse = ", "), ")); map in the ",
    "geom instead, as geom_jbrowse_reads(aes(fill = mapq))"),
    .frequency = "regularly", .frequency_id = "ggjbrowse_plot_mapping")
}

#' @exportS3Method ggplot2::fortify
fortify.jbrowse_result <- function(model, data, ...) {
  placed <- on_axis(model)
  structure(placed[[attr(model, "primary")]], jbrowse = placed)
}

# A ggplot scale built from its defaults, each argument the caller names
# taking the default's place
scale_with <- function(scale, defaults, ...) {
  own <- list(...)
  defaults[names(own)] <- own
  do.call(scale, defaults)
}

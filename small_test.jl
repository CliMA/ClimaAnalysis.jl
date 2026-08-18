time = Float64.(collect(1:(12 * 10)))
lat = Float64.(collect(-90:90))
lon = Float64.(collect(-180:180))
var =
    TemplateVar() |>
    add_dim("time", time, units = "s") |>
    add_dim("lat", lat, units = "degrees") |>
    add_dim("lon", lon, units = "degrees") |>
    add_attribs(start_date = Dates.DateTime(2010), short_name = "lwu") |>
    initialize

module ClimaExplorer

import ClimaAnalysis

import Bonito:
    Observable,
    App,
    Slider,
    Dropdown,
    DOM,
    Asset,
    TextField,
    Card,
    Grid,
    Label,
    Styles,
    Centered,
    Button,
    on
import WGLMakie
import GeoMakie

include("layouts.jl")

function _logo_img()

    # TODO: The logo does not resize correctly
    imstyle = Styles(
        "display" => :block,
        "position" => :relative,
        "width" => "180px",
        "max-width" => :none, # needs to be set so it's not overwritten by others
    )
    img = DOM.a(
        href = "https://www.github.com/CliMA/ClimaAnalysis.jl",
        DOM.img(;
            src = Asset(joinpath("assets", "logo-full.svg")),
            style = imstyle,
        ),
    )
    return img
end



"""
    BonitoApp(path)

Return a `Bonito` `App` for data at `path`.
"""
function BonitoApp(path)
    app = App(; title = "ClimaExplorer") do

        # Text field to change the path
        #
        # TODO: width = 98% is pretty random. We should find a truly responsive with.
        path_obs = TextField(path, style = Styles("width" => "98%"))

        # Most of the logic is implemented in the layout.jl file. Every time path_obs
        # change, the body of the web page is redrawn. This allows us to handle invalid
        # states too.
        body_layout = map(layout, path_obs)

        img = _logo_img()
        input_field = Centered(
            Grid(
                "Insert the path of your output data and press ENTER",
                path_obs,
            ),
        )

        grid = Grid(Card(input_field), Card(img); columns = "90% 10%")

        return Grid(grid, body_layout)
    end
    return app
end

end

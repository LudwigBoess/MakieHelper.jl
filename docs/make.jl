using Documenter
using MakieHelper

makedocs(
    sitename = "MakieHelper",
    format = Documenter.HTML(),
    modules = [MakieHelper]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#

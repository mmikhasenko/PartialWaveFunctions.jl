module PartialWaveFunctions

export f_logfact
include("factorials.jl")

export wignerd, wignerd_doublearg
export wignerD, wignerD_doublearg
export kronecker
include("wignerd.jl")

export CG, CG_doublearg
export clebschgordan, clebschgordan_doublearg
include("clebsch_gordan.jl")

end # module

# # # # # # # # # # # # # # # # # # # # # # # #
# created with PkgTemplates
# using PkgTemplates
# t = Template(;
#     user="mmikhasenko",
#     license="MIT",
#     authors="Misha Mikhasenko",
#     dir=joinpath(DEPOT_PATH[1], "dev"),
#     julia_version=v"1.2",
#     plugins=[
#         TravisCI(),
#         Codecov(),
#         AppVeyor(),
#         ],
#     )
# generate("PartialWaveFunctions", t)
# # # # # # # # # # # # # # # # # # # # # # # # # #

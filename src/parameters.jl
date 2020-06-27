Base.@kwdef mutable struct Parameters
    compressor
end

const DEFAULT_PARAMETERS = Parameters(compressor = PartialACA(atol=1e-5)
                                      )

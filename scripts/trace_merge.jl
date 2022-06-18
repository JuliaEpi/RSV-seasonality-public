using ArgParse
using DrWatson
include(scriptsdir("utils.jl"))


function trace_merge(runs, ignorecommit=false)

    # Load the arguments early so we can checkout the correct version.
    runs_args = map(runs) do rundir
        DrWatson.wload(joinpath(rundir, "args.jld2"))
    end;

    # Checks
    @info "(♻) Comparing model runs for consistency..."

    if length(unique((x["betafunc"] for x in runs_args))) > 1
        println("(⚠) Oh no!!! betafunc is not identical among all runs. Exit")
        exit(1)
    end

    if length(unique((x["rhofunc"] for x in runs_args))) > 1
        println("(⚠) Oh no!!! rhofunc is not identical among all runs. Exit")
        exit(1)
    end


    if length(unique((x["differential-equation"] for x in runs_args))) > 1
        println("(⚠) Oh no!!! the ODE function is not identical among all runs. Exit")
        exit(1)
    end

    unique_gitcommit = unique((x["gitcommit"] for x in runs_args))
    if length(unique_gitcommit) > 1
        println("(⚠) Oh no!!! Runs come from different commits; do you want to continue? [y/N]: ")
        answer = readline()
        if lowercase(answer) != "y"
            println("Exiting...")
            exit(1)
        end

        gitcommit = runs_args[1]["gitcommit"]
        println("Okay then. We'll check out $(gitcommit) then.")
    end

    if !ignorecommit
        unique_gitpatch = unique((getkey(x, "gitpatch", "") for x in runs_args))
        if length(unique_gitpatch) > 1
            @warn "(⚠) Oh no!!! Runs come come from different patches. Applying the first one..."
        end

        @info "(♻) Checking out $(first(unique_gitcommit))..."
        checkout!(runs_args[1])
        @info "(✓) Checkout completed!"
    end

    # Load the traces.
    traces = map(runs) do rundir
        deserialize(joinpath(rundir, "trace.jls"))
    end;

    # Load the models.
    models = map(runs) do rundir
        deserialize(joinpath(rundir, "model.jls"))
    end;

    # Combine the chains if so specified.
    @warn "Assuming first args/model is representative of all args/models."
    run_args = runs_args[1]
    model = models[1]

    @info "Combining $(length(traces)) chains into a single chain."
    trace = MCMCChains.chainscat(traces...)

    # Drop the warmup.
    nadapts = run_args["nadapts"] 
    trace = trace[nadapts + 1:end] 
    @info "Dropping $(nadapts) adaptation steps; $(length(trace)) samples remaining."

    return model, trace, run_args
end


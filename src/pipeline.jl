module Pipeline

# Define an abstract type since I want define one method to get what variable the pipeline
# is defined for

# Define a marco that make a struct with a specified name (e.g. ObsPipeline or SimPipeline)
# and subtype the abstract type
# I think this should just contain a field of a vector of short name
# Should also define a generic method that will fail on any OutputVar whose pipeline
# is not defined

# Idea is that we can form a Val type from the short name of the OutputVar, so we can do
# stuff like pipeline(var) and it will pick up the right pipeline to operate on

# Define a marco that is of the form
# @pipeline "short_name1" "short_name2" ... arg1 arg2 arg3 kwarg1 kwarg2 begin
# body
# end
# One thing that might be annoying is having to define the arguments and kwargs every
# single time, but maybe it should be like that since you might want to modify a part of the pipeline without knowing what the other stuff do?

# Also be able to define a pipeline with no short name, so that it is a generic pipeline that
# can take in OutputVar with no name too!



end

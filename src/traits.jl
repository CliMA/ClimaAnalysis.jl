"""
    abstract type HasDimensionMetadataTrait end

An object that indicates whether something contains information about the attributes and
dimensions.

For any object `x` where `HasDimensionMetadataTrait(typeof(x)) = true`, the dimensions,
attributes, and dimension attributes must be accessible by `x.attributes`, `x.dims`, and
`x.dim_attributes`. These fields must be implemented identically to `OutputVar`.
"""
abstract type HasDimensionMetadataTrait end

"""
    HasDimMetadata

The object's attributes, dimensions, and dimension attributes can be accessed through
fields `attributes`, `dims`, and `dim_attribs`.
"""
struct HasDimMetadata <: HasDimensionMetadataTrait end

"""
    HasNoDimMetadata

The object does not contain attributes, dimensions, or dimension attributes or it cannot be
accessed through the fields `attributes`, `dims`, and `dim_attribs`.
"""
struct HasNoDimMetadata <: HasDimensionMetadataTrait end

# Define trait for only OutputVar, FlatVar, and Metadata.
# This let us reuse functions that does not require accessing the data. For example, you may
# want to access the longitudes, latitudes, or units of a `FlatVar` or `Metadata`.
#
# An alternative design is to define a global constant that is an union of all the types.
#
# TODO: Ask Nat about this since I am not sure if this is overkill or if it is easier to
# just use an union of types (e.g. const HasDimMetadata = Union{OutputVar, FlatVar, Metadata})
# The difference is extensibility from the perspective of the user, since you can easily
# modify the const union of types
HasDimensionMetadataTrait(::Type{<:OutputVar}) = HasDimMetadata()
HasDimensionMetadataTrait(::Type{<:FlatVar}) = HasDimMetadata()
HasDimensionMetadataTrait(::Type{<:Metadata}) = HasDimMetadata()
HasDimensionMetadataTrait(::Type) = HasNoDimMetadata()

# Example
# longitudes(x::T) where {T} = longitudes(HasDimensionMetadataTrait(T), x)
# longitudes(::HasDimMetadata, x) = [0.0, 1.0, 2.0]

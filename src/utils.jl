## Common utility functions

# function _custom_metafmt(level::Logging.LogLevel, _module, group, id, file, line)
#     color = Logging.default_logcolor(level)
#     prefix = (level == Warn ? "Warning" : string(level))*':'
#     return color, prefix, ""
# end
# Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Info, _custom_metafmt, true, 0, Dict{Any,Int}()))

# function ifwarn(msg, ::Val{T}) where T
#     global DOWARN
#     if DOWARN[]
#         with_logger(customLogger::Logging.ConsoleLogger) do
#            @logmsg T msg
#         end
#     end
#     nothing
# end
# ifwarn(msg) = ifwarn(msg, Val(Warn))

macro ifwarn(ex)
    return quote
        if (DOWARN[]::Bool)
            $(esc(ex))
        end
    end
end

function ifexport(c, _name=nothing, path=tempdir())
    global DOEXPORT
    if (DOEXPORT[]::Bool)
        name::String = _name isa Nothing ? tempname(path; cleanup=false)*".vtf" : begin
            i = 0
            x = "crystal_"*_name*'_'*string(i)*".vtf"
            names = readdir(path; join=false, sort=true)
            j = searchsortedfirst(names, x)
            while j <= length(names) && names[j] == x
                i += 1
                j += 1
                x = "crystal_"*_name*'_'*string(i)*".vtf"
            end
            x
        end
        truepath = joinpath(path, name)
        @info "Automatic export of input is enabled: saving file at $truepath"
        export_vtf(truepath, c, 6)
    end
    nothing
end

"""
    soft_widen(::Type)

Internal function used to selectively widen small integer and rational types.
This is useful to avoid overflow without sacrificing too much efficiency by
always having to resolve to very large types.
"""
soft_widen(::Type{T}) where {T} = T
soft_widen(::Type{Int32}) = Int64
soft_widen(::Type{Int16}) = Int32
soft_widen(::Type{Int8}) = Int16
soft_widen(::Type{Rational{T}}) where {T} = Rational{soft_widen(T)}

"""
    issingular(x::SMatrix{3,3,T,9}) where T<:Rational

Test whether a 3x3 matrix is singular. The input matrix must have a wide enough
type to avoid avorflows. If an overflow happens however, it is detected and results
in an error (rather than silently corrupting the result).
"""
function issingular(x::SMatrix{3,3,T,9})::Bool where T<:Rational
    (i, j, k) = iszero(x[1,1]) ? (iszero(x[1,2]) ? (3,1,2)  : (2,1,3)) : (1,2,3)
    iszero(x[1,i]) && return true
    U = widen(T)
    x1i = U(x[1,i])
    factj = (x[1,j] // x1i)
    factk = (x[1,k] // x1i)
    y11 = x[2,j] - factj * x[2,i]
    y12 = x[3,j] - factj * x[3,i]
    y21 = x[2,k] - factk * x[2,i]
    y22 = x[3,k] - factk * x[3,i]
    try
        return y11 * y22 == y12 * y21
    catch e
        @assert e isa OverflowError
        return widemul(y11, y22) == widemul(y12, y21)
    end
    # This can overflow so the input matrix should already have a wide enough type
end

function isrank3(x::AbstractMatrix{T}) where T<:Rational
    _n, n = size(x)
    @assert _n == 3
    n < 3 && return false
    cols = collect(eachcol(x))
    u1 = popfirst!(cols)
    u2 = u1
    while iszero(u1) && !isempty(cols)
        u1 = popfirst!(cols)
    end
    n = length(cols)
    k = iszero(u1[1]) ? iszero(u1[2]) ? 3 : 2 : 1
    i = 1
    while i < n
        u2 = cols[i]
        widen(u2[k]//u1[k])*u1 == u2 || break
        i+=1
    end
    i == n && return false
    for j in i+1:n
        a = SMatrix{3,3,T,9}(hcat(u1, u2, cols[j]))
        issingular(a) || return true
    end
    return false
end

"""
    back_to_unit(r::Rational)

Return `x - floor(Int, x)`
"""
function back_to_unit(r::Rational)
    return Base.unsafe_rational(mod(numerator(r), denominator(r)), denominator(r))
end


"""
    nextword(l, k)

Return the triplet of indices `(i, j, x)` such that `l[i:j]` is the next word in the
string `l` after position `k`.
Use `k = x` to get the following word, and so forth.

`(0, 0, 0)` is returned if there is no next word.
"""
function nextword(l, i)
    n = lastindex(l)
    i == n && return (0, 0, 0)

    i::Int = nextind(l, i)
    start::Int = 0
    while i <= n
        if isspace(l[i])
            i = nextind(l, i)
        else
            if l[i] == '#'
                i = findnext(isequal('\n'), l, i)
                i == nothing && return (0, 0, 0)
            else
                start = i
                break
            end
        end
    end
    start == 0 && return (0, 0, 0)

    inquote = l[start] == '\'' || l[start] == '"'
    quotesymb = l[start]
    inmultiline = l[start] == ';' && l[prevind(l,start)] == '\n'
    instatus = inmultiline | inquote
    instatus && (start = nextind(l, start); i = start)
    while i <= n
        c = l[i]
        if !isspace(c)
            if instatus
                if inmultiline && c == ';' && l[prevind(l,i)] == '\n'
                    return (start, prevind(l, i-1), i)
                elseif inquote && c == quotesymb && i != n &&
                        (isspace(l[nextind(l, i)]) || l[nextind(l, i)] == '#')
                    return (start, prevind(l, i), i)
                end
            elseif c == '#'
                return (start, prevind(l, i), prevind(l, i))
            end
            i = nextind(l, i)
        elseif !instatus
            return (start, prevind(l, i), i)
        else
            i = nextind(l, i)
        end
    end
    if instatus
        if inmultiline
            throw("Invalid syntax: opening multiline field at position $start is not closed")
        end
        if inquote
            throw("Invalid syntax: opening quote $quotesymb at position $start is not closed")
        end
    end
    return (start, n, n)
end

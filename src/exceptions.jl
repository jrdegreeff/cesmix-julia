# Generally I'm not a fan of custom exceptions. It adds quite some bloat
# and in this case, why would that be better than getting a Julia stacktrace?
#
# Exceptions

struct UnimplementedError <: Exception
    func::Symbol
    type::Type
    UnimplementedError(func::Symbol, any) = new(func, typeof(any))
end

Base.showerror(io::IO, e::UnimplementedError) = print(io, e.type, " doesn't implement ", e.func)

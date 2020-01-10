import Base: *
import Base: +
import Base: -
import Base: /
import Base: ^
import Printf.@printf
import Printf.@sprintf
using Calculus

struct Error
    A
    a
end

Error(x) = Error(x,0)

degToRad(x) = 2*π*(x/360)
radToDeg(x) = 360*(x/(2*π))
I(x) = x

function *(x, y::Error)
    Anew = x*y.A
    anew = x*y.a
    @printf("%.8g(%.8g \\pm %.8g) = %.8g \\pm %.8g\n\n", x, y.A, y.a, Anew, anew)
    return Error(Anew, anew)
end

function +(x::Error, y::Error)
    Anew = x.A + y.A
    anew = (x.a^2 + y.a^2)^0.5
    @printf("(%.8g \\pm %.8g) + (%.8g \\pm %.8g) = (%.8g + %.8g) \\pm \\sqrt{%.8g^2 + %.8g^2} = %.8g \\pm %.8g\n\n", x.A, x.a, y.A, y.a, x.A, y.A, x.a, y.a, Anew, anew)
    return Error(Anew, anew)
end

function -(x::Error, y::Error)
    Anew = x.A - y.A
    anew = (x.a^2 + y.a^2)^0.5
    @printf("(%.8g \\pm %.8g) - (%.8g \\pm %.8g) = (%.8g - %.8g) \\pm \\sqrt{%.8g^2 + %.8g^2} = %.8g \\pm %.8g\n\n", x.A, x.a, y.A, y.a, x.A, y.A, x.a, y.a, Anew, anew)
    return Error(Anew, anew)
end

function *(x::Error, y::Error)
    Anew = x.A * y.A
    anew = Anew * ((x.a/x.A)^2 + (y.a/y.A)^2)^0.5
    @printf("(%.8g \\pm %.8g)(%.8g \\pm %.8g) = (%.8g \\cdot %.8g) \\pm %.8g \\cdot %.8g\\sqrt{\\left(\\frac{%.8g}{%.8g}\\right)^2 + \\left(\\frac{%.8g}{%.8g}\\right)^2} = %.8g \\pm %.8g\n\n", x.A, x.a, y.A, y.a, x.A, y.A, x.A, y.A, x.a, x.A, y.a, y.A, Anew, anew)
    return Error(Anew, anew)
end

function /(x::Error, y::Error)
    Anew = x.A / y.A
    anew = Anew * ((x.a/x.A)^2 + (y.a/y.A)^2)^0.5
    @printf("\\frac{%.8g \\pm %.8g}{%.8g \\pm %.8g} = \\frac{%.8g}{%.8g} \\pm \\frac{%.8g}{%.8g}\\sqrt{\\left(\\frac{%.8g}{%.8g}\\right)^2 + \\left(\\frac{%.8g}{%.8g}\\right)} = %.8g \\pm %.8g\n\n", x.A, x.a, y.A, y.a, x.A, y.A, x.A, y.A, x.a, x.A, y.a, y.A, Anew, anew)
    return Error(Anew, anew)
end

function ^(x::Error, y)
    Anew = x.A ^ y
    anew = x.A^(y-1) * y * x.a
    @printf("(%.8g \\pm %.8g)^%.8g = %.8g^%.8g \\pm %.8g^%.8g \\cdot %.8g \\cdot %.8g = %.8g \\pm %.8g\n\n", x.A, x.a, y, x.A, y, x.A, y-1, y, x.a, Anew, anew)
    return Error(Anew, anew)
end

function funEval(f::Function, g::Function, x::Error, useRaw::Bool = false)
    Anew = f(g(x.A))
    anew = derivative(f, g(x.A)) * x.a
    if(useRaw)
        funName = string(typeof(f).name.mt.name)
        @printf("\\%s(%.8g \\pm %.8g) = \\%s(%.8g) \\pm %.8g\\cdot\\frac{df}{dA} = %.8g \\pm %.8g\n\n", funName, x.A, x.a, funName, x.A, x.a, Anew, anew)
    end
    if(!useRaw)
        funName = "f"
        @printf("%s(%.8g \\pm %.8g) = %s(%.8g) \\pm %.8g\\cdot\\frac{df}{dA} = %.8g \\pm %.8g\n\n", funName, x.A, x.a, funName, x.A, x.a, Anew, anew)
    end
    return Error(Anew, anew)
end

function combineData(data::Array{Float64}, varName::String)
    ave = 0
    for i in data
        ave += i
    end
    ave /= length(data)

    @printf("\\bar{%s} = \\frac{1}{%.8g}\\sum_{i = 1}^{%.8g}%s_i = %.8g\n\n", varName, length(data), length(data), varName, ave)
    varr = zeros(Float64, length(data))

    cancer1 = ""
    for i in length(data)
        newStr = @sprintf("(%.8g - %.8g)^2", ave, data[i])
        if i < length(data)
            cancer1 = string(cancer1, newStr, " + ")
        end
        if i == length(data)
            cancer1 = string(cancer1, newStr)
        end
    end

    cancer2 = ""
    for i in length(data)
        varr[i] = (ave - data[i])^2
        newStr = @sprintf("(%.8g)^2", ave - data[i])
        if i < length(data)
            cancer2 = string(cancer1, newStr, " + ")
        end
        if i == length(data)
            cancer2 = string(cancer1, newStr)
        end
    end

    fullVarr = 0;
    for i in length(data)
        fullVarr += varr[i]
    end

    stdev = (fullVarr/(length(data) - 1))^0.5

    error = stdev/(length(data)^0.5)

    finalRes = string("\\sigma = ", @sprintf("\\sqrt{\\frac{\\sum_{i=1}^{%.8g}(\\bar{%s} - %s_i)^2}{N - 1}}", length(data), varName, varName), " = ", @sprintf("\\sqrt{\\frac{%s}{%.8d}}", cancer1, length(data)-1), " = ", @sprintf("\\sqrt{\\frac{%s}{%.8d}}", cancer2, length(data)-1), " = ", @sprintf("\\sqrt{\\frac{%.8d}{%.8d}}", fullVarr, length(data)-1), " = ", stdev)
    println(finalRes)
    println()

    lastOne = string("\\frac{\\sigma}{\\sqrt{N}} = ", error, "\n")
    println(lastOne)

    final = @sprintf("%.8d \\pm %.8d\n", ave, error)
    println(final)

    return Error(ave, error)
end

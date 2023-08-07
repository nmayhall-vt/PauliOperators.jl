struct MyType{T}
    a::T
    b::T
end


function Base.:*(mt1::MyType{T}, mt2::MyType{T}) where T
    return MyType{T}(mt1.a * mt2.a, mt1.b * mt2.b)
end


function test(T)
    a = MyType{T}(1,2)
    b = MyType{T}(3,4)
    @btime c = $a * $b
end

function testTuple()
    a = (1,2)
    b = (3,4)
    @btime c = ($a[1] * $b[1], $a[2] * $b[2])
end
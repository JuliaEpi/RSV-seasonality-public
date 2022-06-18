function NegativeBinomial2(ψ, increp)
    p = 1.0/(1.0 + ψ*increp)
    #p = 1.0/(1.0 + increp/ψ)
    r = 1.0/ψ
    #r = ψ
    return NegativeBinomial(r, p)
end

function nbinomlogpdf2(ψ, increp, x)
    p = 1.0/(1.0 + ψ*increp)
    r = 1.0/ψ

    return StatsFuns.nbinomlogpdf(r, p, x)
end

function uiyx(ui, p)
    N = length(p[:, 1])
    uix = ui .* p[Int(N/2), :]
    uiyx = repeat(uix', N, 1)
    return uiyx
end
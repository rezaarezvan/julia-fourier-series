# Archimedian spirals
using Luxor

path   = "C:/Users/Reza/Desktop/juliatesting/media/archimedian.png" # Absolute, change this
heigth = 600
width  = 600

spiraldata = [
  (-2, "Lituus",      130),
  (-1, "Hyperbolic", 180),
  ( 1, "Archimedes",   1),
  ( 2, "Fermat",       7)]

grid = GridRect(O - (200, 0), 130, 50)

@png begin
    background("grey7")
    for aspiral in spiraldata
        @layer begin
            sethue("white")
            setcolor(sethue("white"))
            setline(1)
            Luxor.translate(nextgridpoint(grid))
            spiral(last(aspiral), first(aspiral), period=15π, :stroke)
            fontsize(20)
            label(aspiral[2], :N, offset=80)
        end
    end
end heigth width path 



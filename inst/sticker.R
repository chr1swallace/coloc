## devtools::load_all("~/RP/coloc")
library(colorspace)
library(data.table)
library(grid)
library(coloc)
data(coloc_test_data)
attach(coloc_test_data)
res=coloc.abf(D3,D4)
library(ggplot2); library(cowplot); theme_set(theme_cowplot(font_size=10))
library(magrittr)
library(seaborn)
show_seaborn()

set.seed(42)
df=data.table(x=res$results$position,
              y1=-log10(pnorm(-abs(res$results$z.df1))*2),
              y2=ifelse(res$results$position < 50,
                        -log10(pnorm(-abs(res$results$z.df1))*2),
                        -log10(pnorm(-abs(res$results$z.df2))*2)))
df[,mag:=x < 50 & y1 > 5]

col1=seaborn:::SEABORN_PALETTES$muted6[c(5)] # background and border
col2=seaborn:::SEABORN_PALETTES$deep6[c(6)] # text and highlight points

col1="#4831d4"
col2="#ccf381"

## orange/wheat
col2=darken("#ee4e34",0.2)
## col1="#FBF3EA" #"#fcedda"
col1="#fcedda"

## ## yellow/blue
## col2="#234e70"
## col1="#fbf8be"

## ## yellow/blue/green
## col1="#ffd55a"
## col2="#293250"
## col2="#66d45e"

set.seed(42)
p=ggplot() +
  geom_point(aes(x=x+3,y=y1+1),col="grey70",data=df[mag==FALSE],size=1) +
  geom_jitter(aes(x=x,y=y2),col="grey30",data=df[mag==FALSE],size=1) +
  geom_point(aes(x=x+3,y=y1+1),col=lighten(col2,0.4),data=df[mag==TRUE],size=1.5) +
  geom_jitter(aes(x=x,y=y2),col=darken(col2,0),data=df[mag==TRUE],size=1.5) +
  ylim(0,12) +
  background_grid(major="y") +
  theme(axis.text=element_blank(),axis.title=element_blank(),axis.line.y=element_blank())
  ## p

x=c(0.85,0.455)
y=c(0.05,0.525)
w=0.9
xm=(1-w)*x[1] + w*x[2]
ym=(1-w)*y[1] + w*y[2]
w2=.35
xmag=x[2] + w2*diff(x)
ymag=y[2] + w2*diff(y)

circo=circleGrob(x=xmag, y=ymag, r=0.25, default.units="npc", name=NULL,
           gp=gpar(lwd=3,col=darken(col2,0.0)), vp=NULL)
lineo=linesGrob(x = unit(c(x[1], xm), "npc"),
          y = unit(c(y[1], ym), "npc"),
          default.units = "npc",
          arrow = NULL, name = NULL,
          gp=gpar(lwd=4,col=darken(col2,0.0)), vp = NULL)
linei=linesGrob(x = unit(c(xm, x[2]), "npc"),
          y = unit(c(ym, y[2]), "npc"),
          default.units = "npc",
          arrow = NULL, name = NULL,
          gp=gpar(lwd=2,col=darken(col2,0.0),lineend="butt"), vp = NULL)
p2=ggdraw(p) +
  ## draw_grob(circo,scale=0.5,x=-0.2,y=0.2) +
  ## draw_grob(circo,scale=0.48,x=-0.2,y=0.2) +
  draw_grob(circo) +
  ## draw_grob(circo,scale=0.96) +
  draw_grob(lineo) +
  draw_grob(linei)
## p2

library(hexSticker)
s <- sticker(p2,
             package="coloc",
             p_size=24, p_color=col2, p_fontface="bold", #darken(col2,0.2),
             p_y=1.5,
             s_x=1, s_y=.8, s_width=1.4, s_height=1,
             u_size = 1,
             h_fill=darken(col1,0.00),#lighten(col1,0),
             h_color=darken(col1,0.2),
             filename="~/RP/coloc/man/figures/logo.png")
plot(s)
system("magick convert ~/RP/coloc/man/figures/logo.png -resize 50% ~/RP/coloc/man/figures/logo50.png")
system("magick convert ~/RP/coloc/man/figures/logo.png -resize 30% ~/RP/coloc/man/figures/logo30.png")

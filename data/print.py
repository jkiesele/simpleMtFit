

numbers=[
        830.91  ,      
    830.167     ,   
    831.758     ,   
    831.279     ,   
    832.913     ,   
    830.346     ,   
    828.873     ,   
    829.571     ,   
    828.985     ,   
    832.387     ,   
     830.281    ,    
     829.257    ,    
     831.195    ,    
     830.051    ,    
     831.728    ,    
     831.418    ,    
     831.884    ,    
     830.378    ,    
     830.782    ,    
     831.144    ,    
     831.57     ,   
     831.086    ,    
     830.943    ,    
     831.831    ,    
     829.734    ,    
     831.423    ,    
     832.209    ,    
     829.646    ,    
     829.556    ,    
     834.908    ,    
     831.183    ,    
     831.568    ,    
     831.435    ,    
     831.456    ,    
     831.579    ,    
     830.865    ,    
     830.724    ,    
     831.497    ,    
     831.033    ,    
     830.452    ,    
     831.102    ,    
     830.654    ,    
     829.636    ,    
     830.962    ,    
     831.648    ,    
     830.57     ,   
     831.883    ,    
     831.753    ,    
     830.043    ,    
     830.57     ,   
     831.227    ,    
     831.436    ,    
     830.605    ,    
     831.071    ,    
     831.207    ,    
     831        ,
     830.923    ,    
     830.961    ,    
     831.074    ,    
     830.687    ,    
     830.866    ,    
     830.6      ,  
     830.803    ,    
     830.927    ,    
     830.734    ,    
     831.04     ,   
     831.062    ,    
     830.819    ,    
     831.153    ,    
     831.122    ,    
     830.264    ,    
     830.428    ,    
     831.011    ,    
     830.516    ,    
     830.479    ,    
     830.538    ,    
     831.258    ,    
     830.03     ,   
     831.094    ,    
     831.159    ,    
     831.216    ,    
     831.137    ,    
     831.247    ,    
     831.005    ,    
     830.732    ,    
     831.13     ,   
     831.096    ,    
     830.484    ,    
     831.038    ,    
     830.77     ,   
     831.055    ,    
     830.753    ,    
     830.828    ,    
     831.173    ,    
     830.845    ,    
     831.233    ,    
     830.582    ,    
     831.583    ,    
     830.65     ,   
     830.924    ,    
      830.891  
]


#no resumm
numbers=[
       811.89  ,
   811.167     ,
   812.697     ,
   812.228     ,
   813.858     ,
   811.352     ,
   809.901     ,
   810.6       ,
   810.004     ,
   813.344     ,
    811.273    ,
    810.225    ,
    812.135    ,
    810.982    ,
    812.675    ,
    812.375    ,
    812.869    ,
    811.364    ,
    811.733    ,
    812.124    ,
    812.515    ,
    812.08     ,
    811.911    ,
    812.783    ,
    810.727    ,
    812.385    ,
    813.159    ,
    810.651    ,
    810.563    ,
    815.821    ,
    812.151    ,
    812.537    ,
    812.411    ,
    812.429    ,
    812.544    ,
    811.852    ,
    811.706    ,
    812.463    ,
    812.002    ,
    811.442    ,
    812.08     ,
    811.64     ,
    810.637    ,
    811.925    ,
    812.615    ,
    811.556    ,
    812.845    ,
    812.715    ,
    811.037    ,
    811.555    ,
    812.199    ,
    812.392    ,
    811.577    ,
    812.047    ,
    812.183    ,
    811.979    ,
    811.903    ,
    811.925    ,
    812.052    ,
    811.671    ,
    811.845    ,
    811.586    ,
    811.786    ,
    811.906    ,
    811.717    ,
    812.016    ,
    812.029    ,
    811.8      ,
    812.128    ,
    812.1      ,
    811.252    ,
    811.417    ,
    811.99     ,
    811.5      ,
    811.465    ,
    811.522    ,
    812.235    ,
    811.02     ,
    812.071    ,
    812.135    ,
    812.193    ,
    812.113    ,
    812.222    ,
    811.984    ,
    811.714    ,
    812.108    ,
    812.075    ,
    811.469    ,
    812.016    ,
    811.75     ,
    812.034    ,
    811.736    ,
    811.811    ,
    812.151    ,
    811.826    ,
    812.209    ,
    811.565    ,
    812.555    ,
    811.633    ,
    811.905    ,
     811.87
    
    
    ]


lst=[]
sqsum=0
for n in numbers:
    sqsum+=(n-numbers[0])*(n-numbers[0])
    if n-numbers[0] >0:
        print str(n-numbers[0]),
        
print "negative"
for n in numbers:
    if n-numbers[0] <0:
        print str(-(n-numbers[0])),

from math import sqrt

print "sqsum", sqrt(sqsum)
print "rel ", sqrt(sqsum)/numbers[0]

import sys; sys.path.insert(0,'.')
from relabel_overlay import apply, D_mhalf
L='LS'; B='LSB'
# Figure 3 (scope audit): panel D subtitle and legend
apply('F3_orig.pdf','Figure3.pdf',[
 ('white',465,387,621,402),
 ('text',[('Decoy- and entrapment-support indices',L,1,0)],11.2,466.9,399.0,0,'left'),
 ('white',543,406,632,418), ('text',D_mhalf('D'),10,545,415.8,0,'left'),
 ('white',543,418.5,632,430.5), ('text',D_mhalf('E'),10,545,427.8,0,'left'),
])
# Figure 4 (calibration and stability): panel A title, subtitle, right axis; panel D y-axis
apply('F4_orig.pdf','Figure4.pdf',[
 ('white',126,6,215,22), ('text',[('& support index',B,1,0)],13.3,127.1,19.0,0,'left'),
 ('white',122,24.5,230,39), ('text',[('and ',L,1,0)]+D_mhalf('D')+[(' (right axis)',L,1,0)],11.2,123.9,35.6,0,'left'),
 ('white',325.5,108,340.5,140), ('text',D_mhalf('D'),11.2,336.6,124.7,90,'center'),
 ('white',38,305,53,430), ('text',[('D',L,1,0),('b',L,0.75,-0.25),(' / T',L,1,0),('b',L,0.75,-0.25),(' (decoy-derived reference)',L,1,0)],10.5,49.4,368.6,90,'center'),
])
# Figure 5 (MS2Rescore): panel B title and y-axis
apply('F5_orig.pdf','Figure5.pdf',[
 ('white',301,6,467,23), ('text',[('Decoy-support index ',B,1,0),('D',B,1,0),('α','DJBI',0.75,-0.25),('−1/2',B,0.7,0.45),(' vs ',B,1,0),('α','DJBI',1,0)],13.3,302.8,19.2,0,'left'),
 ('white',280,92,296,123), ('text',D_mhalf('D'),11.2,292.0,107.6,90,'center'),
])

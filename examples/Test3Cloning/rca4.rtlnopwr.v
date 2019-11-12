module rca4 (a, b, ci, s, co);

input ci;
output co;
input [3:0] a;
input [3:0] b;
output [3:0] s;

wire vdd = 1'b1;
wire gnd = 1'b0;

BUFX2 BUFX2_1 ( .A(fa0_s), .Y(s[0]) );
BUFX2 BUFX2_2 ( .A(fa1_s), .Y(s[1]) );
BUFX2 BUFX2_3 ( .A(fa2_s), .Y(s[2]) );
BUFX2 BUFX2_4 ( .A(_4_), .Y(s[3]) );
BUFX2 BUFX2_5 ( .A(_0_), .Y(co) );
INVX1 INVX1_1 ( .A(e), .Y(_4_) );
OR2X2 OR2X2_1 ( .A(ci), .B(a[0]), .Y(_5_) );
NAND2X1 NAND2X1_1 ( .A(ci), .B(a[0]), .Y(_6_) );
NAND3X1 NAND3X1_1 ( .A(_4_), .B(_6_), .C(_5_), .Y(_7_) );
NOR2X1 NOR2X1_1 ( .A(ci), .B(a[0]), .Y(_1_) );
AND2X2 AND2X2_1 ( .A(ci), .B(a[0]), .Y(_2_) );
OAI21X1 OAI21X1_1 ( .A(_1_), .B(_2_), .C(b[0]), .Y(_3_) );
NAND2X1 NAND2X1_2 ( .A(_3_), .B(_7_), .Y(fa0_s) );
OAI21X1 OAI21X1_2 ( .A(_4_), .B(_1_), .C(_6_), .Y(c1) );
INVX1 INVX1_2 ( .A(_4_), .Y(k) );
endmodule

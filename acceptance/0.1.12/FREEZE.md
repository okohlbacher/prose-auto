# FREEZE: prose-auto 0.1.12 held-out benchmark (PREREG §10)

- prose-auto commit: bccca5d495b3afa525ce631428792a5c1e84f601
- Jev model: jev-1.14.0
- Jev endpoint: https://api.typesafe.ai/v1/systemone
- Jev votes: 3
- BinAC2 container image (SHA-256, measured on BinAC2): bb7b22210f8588ecb2779d7810ca91beec158ef4446897bc615d9ad4a5a29fef
- pip freeze: /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/pip_freeze.txt

## Step 2: scoring, scenario, draw, power and harness code, manifests, snapshots (2026-11-22 14:14 CET)

```
db23f41e2fd32d0bac5096a81578e5105c2eb112d0cf3b5f4c2977999d3dd8c3  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/PREREG.md
7d3b48a09622ccde90a11312a390489e018809a293f9cbee2b8d3eb75594f485  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/SHA256SUMS
79d72a2294c84a0e73e860da8f8faac59f7dbc2101912474816d7832f82d8a1c  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/extract.py
b78749c4115de33a97918a48dc250a22d104f89c1ebbe9859aa95d19fbd4b2d7  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/score.py
52d846492f4eb3b634e8e98cb7d1b9af656fc14ea74f66d82b48a43a79708193  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/scenarios.py
154d3e7140b601b021b807769d603e7b149f854f2b5faf73883924b340befc55  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/resolve_uniprot.py
8ab10ef6d5c0235f4e7a1cadc6db14b17477f175d1b09b325add9c5a495d2551  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/refdb.py
c191e13307682fe299b929f53ea7bf9faf1c078b731684d756b655de520fb0d3  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power.py
41ab6646831a6091b4f3400f4b42d9cf5d42106c22e14f964debada0ec8a5c10  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power.txt
fb888637ee702cd1a000f905e3cd78377431ded4fe72466ef1c3657ad8f39069  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_grid.txt
2c97609731710ceb6192e57710b3b14f17bc75e4fed1e78f0826b2a65465d770  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_decision.py
0ff31ada93f4182b27e4af896a5658df939395481477d80dc11cbbbf53575316  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_decision.txt
fcf935835e0c536771fc837ac950a74085c905090bdc4f4b3cae67d582f3ab32  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/draw_sets.py
dc9b8d8f52b182a1a523522f483cacc5ad8b1fc1dab38b6ac63e5cd4de444143  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/draw_topup.py
d23cbd0607d9165b68eb38f658868208dccce851bf6b0b97cde31356a33d67b5  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/exclude.txt
93d30a52f1a4df5859c8112d863c73ac6ae1da7f2276953b82621064f01b7ad3  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout.tsv
b249f1ece76e45b40b11d197a74d1da9128a70b2aaf4a9a3e20d22d60e36b24e  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/topup.tsv
a2c8435dc76f2ccf2e7451c8d096650e7092ced89578a08ea1a0e1364a6b8a0b  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/wrongdb.tsv
32f096fe64c8b122208d2910b41825e13450800cc48f161f5ad84647d5618e42  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/candidates.tsv
29fc41dd97099a49ea60365909ffaee4d2caed5ac839bc5e4efcf0ae41e82624  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/topup_candidates.tsv
1c57a4cb50b5bfd9868fd2d01b3122dbaa05ed507c6c14d7d6b1f9dec9d13c5d  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/scenarios.tsv
0eedf45283c2f09ac491524f19afef7a6774cc2d75f22becea691f1f80b63a8e  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/refdb_draw.tsv
2df90d911f8b082feb8eecdb95c2af171f23d24c1508c80131011aca5868de8c  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_inputs/dev_checks.tsv
b8f6e1c782fd8e5d80444d95ad6bafc6076c3bbef9fad0d476108d2007c9cec6  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_inputs/related_ranks.tsv
618aaa9faacdaba43818f6ac9be2c26bec48ff5b0103205fccd2b86d114e32e6  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/power_inputs/related_strata.txt
3b692a148e50ace8b5d10bd441c298a97867ebf418bceb9cf251e44cb9ac448b  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot/candidates.tsv
8821546daeda4e683ac9e627ffbb8f5e30d55cbb6125af3820fd072f5d47e6af  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot/names.dmp
b44e3de08cc5d51f8d48cc038694fcabd609a9758d14b8117d4d8a1eb1541b6d  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot/nodes.dmp
3354f4f4873e8881b5a34f84574f72524f70a9487e4ffa5f9e6488686dd5d837  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot/queries.tsv
60e2e969cbc967f6f072ea2d9eb0c632c48572aa5fc9bc194dc0103f582dd60c  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot_heldout/candidates.tsv
46f2c47a079dd9f00cc06e451ed13398149eabdacf3240f5a3cadf47f90cbd48  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot_heldout/names.dmp
f5640624f346b88f08d89d28b156de65d9823dc098f8d5ddc78c7681385e4d25  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot_heldout/nodes.dmp
3a8f462e410428d51dce5409431403de4f3492568a329e394cea2a0060562efb  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/uniprot_heldout/queries.tsv
d3f8510431cd9e7228d57fb13c4e8f94eaf3d3d3a0f4860e9fa555fb0415c08b  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/annotate.py
c0d47d0f7cdb89ce7be8dd23226dba483bc52bfc4a589daa5ace3e5b262c9278  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/arms.py
0185dfcd3d2204df83ff555be0a7e6dd2f1c5b321097af8914a761947b946ecb  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/freeze.py
681f31bb86c74b1c6a9e03d52cfcb6998a91e406d62bf6de0444440368400e75  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/harness.py
d7c29d59f8602a92f52b36c52a8b09583e5c02eccdf722b99207fd673eca8174  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/heldout.slurm
924013c412ace40e20a8e801157ba7d1cd97e041ffce006dadec7b703a86125e  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/heldout/transfer.py
3b714d55680d34285caaa77671c800928393d6cf917d3f4127392284ab3bd5e2  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/corpus/corpus.py
d259efe7fef73e3af4ffe5aee474b1f40cad0b31b9b7155a8cecdaafca080d4e  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/corpus/entrap.py
bfa9dc5bf03c8252fb6caa0e58d5af7062d8ea73f8e6284445b639e3862164a0  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/acceptance/0.1.12/corpus/entrap_pre22c.py
```

## Step 2: FASTA files (UniProt proteomes, the reference panel, cRAP) (2026-11-22 14:14 CET)

```
3d63b778aa8bfeae49d5a8bc5096ce18707229e1a9049341a720ed9d27f2712a  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000437.fasta
dbf259de0c637fb24010152bdf7249af733d3e4830ded9e953d39bc679edd0df  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000539.fasta
5e49f15913722d518ddc7662fb2fd663de6ac72e8200c1dd7c202b70e79694e9  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000559.fasta
bc7b19b38120bcbfdd6ca9548b3f7545dd9602f66f3f0a273d238f9e6f99da62  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000585.fasta
4f3e21581a1b00173ff992f9ebba26a6918f964af1ca67d073551acc96c87c9c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000589.fasta
c67bb322c2d9c963b7e6d6fc260e8a59048a9de56ad67c9a856a4357dba29035  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000625.fasta
65108755b1afec131546a09cfe8ebb6168a143c9e9f4c04497620a6c5bce4b5e  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000653.fasta
751bffc6f9718e75f41b25396b937814229169efa87cd5f47018e9075698708d  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000803.fasta
d98bc3865d6f8a76314ec098774e7deea323484aa36f92390f5191ff86eac88d  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000000817.fasta
83aa9913c0c86de19d7cb3f786bb918fe3249ebc3d2a8a2868e5c01f6c2694ca  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001172.fasta
20fc109d91ed4b834bd26cc1fc04d727b6fb04ec2ed28ef2ffd81b9dc15c0e47  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001410.fasta
03cfbfaa623e4969186bea40d2f10e0a86ca44a1c06b90a5717d5d436838ecc2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001450.fasta
56c85cdd4cca5663f7e439a666769016aaec05bf733f7d189ae67b6218faff5c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001798.fasta
c640d77b8237453d5df0c2a94eccdd10eca76cae7629e2f40c40bd9c250c1519  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001811.fasta
4e1ca0066dd2fb34c12d10f62b4cdd4ff3bce6415994bb9fb053b73c75af0793  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000001940.fasta
f5bbfb0f41cd235f701d6dfc501d7c8a85ba2772d05c6bced23ed2ace52569ed  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000002277.fasta
ff67ac6e6c9192838c1705199ccd9e9c751f372e22264846e9a39e8756e7a105  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000002311.fasta
16882ed6dc4b9d472fd0d034be8790556463a149a241b473701a07076c90d436  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000002485.fasta
7bfaa4e8d832e2bc09bf5eebac9765b7d386e1e9a3ffd846c47448d42b9c5257  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000002494.fasta
ca378270619c48caa9722d83b30f305bf1d7785df57d7a40852a97a036e307ac  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000002564.fasta
8cf52ef95421dca415e03d3d5a6b698a464b0daf8d1caa1368706fa8be43aef3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000004994.fasta
d3b612b4718ed1dc6924b33cc637e89e91ee53c3c3bf757e84ba4169d19907fd  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000005222.fasta
6f649d9224beabad4285cbed762455367f5d0387c0862bacd7314449c517f417  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000005640.fasta
4f33e78803ced5114584d5f1c7a0ce8cdf83c5534596eef5e61a4eb3f7ffc348  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000006548.fasta
e727da43fa35a3f384537cf9324aaaf27654aa3101ce4293fd53839d3c48f387  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000006906.fasta
5083aeb71fc4851f8a17e3e3903d1ac969ee47ebb87c7cb9db84ce7bbdc829bf  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000007305.fasta
8781104df9e8dd1cdf53b89292e855c9af272d48291dfba10747ab6b0cee309c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000008144.fasta
091c88de68f45f65056b9fb128245363f70915bc16e539a08f8ac0adf58c80f8  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000008227.fasta
b9cc75c92fe801f89fef36e5726b64b9d771441b53f69504a0c1f64b5c45381c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000008524.fasta
55f8c677af2d5ffceeae7f7fabbbf0f94d35c417e6bae1467fdb457a7be9f854  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000008816.fasta
90be1ce93705f65517cfe3fbc9789648debe7abf22ee4bc11bf7aad27fdd1ed0  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000009136.fasta
bf2de2111d2ddd741442afa944171ce65a4d5f4aa84c0913dc8376df6781525a  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000011115.fasta
795d1caf535c43d007fdac3630ac5c5897016d9d7e73212913c20522493c11f2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000011116.fasta
1691716754e1226e946200c67e9c21d9514c167b0a74bd133422cda36d4f6f03  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000019116.fasta
a86b0e7c4a583eaf40fa4a82793c28c58e89e8757fd96bc4dd9418784cdcfd6c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000053958.fasta
19b714addca23259b51004c2b663570165e28395ac965ea0b01ed1a7f27b9b40  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000069940.fasta
ba52f264befb6eeaf9741078ed2354dddd761751c61e2e24c68da8c411be1e24  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000186698.fasta
195c62b3039fae76746d1d0556b5f2e9bbc03294511ab3010450b857013dfcea  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000233100.fasta
7f1a2c15843c1c3acfc20164ba7c305cda6aa80429c18b7aef6a349b6da32ceb  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000236500.fasta
d187aa22b8cd239d730867746b038f9ed154f6063a4a552a5472d90cafa48198  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000244005.fasta
ee43e442f799fdfef875ef7c4c8a6b7c10edc9e4b63208c590fa63389f94dac3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000249293.fasta
33ffd9a1b70ff19c0f510a8a6aef3b31c298c5c9b78c097c49c97d6ffb119b35  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000254640.fasta
136ff88e0c05494f6c9cd09a8c1e0129076b6e0802276601a7940c6e9a5e1d56  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000257200.fasta
af263c7cd8ebc5bda1f121ff9b08d748486d0810ee7b1ea2672a50f4c9c216b7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000298196.fasta
0747e88316397f640632153c5a740fbb486074cec72ca4c1b1177f79f664167c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000318052.fasta
95b6bc4a0134a9054351548ef2e0935ee6fd45b54b5f9ec41856c6802923eba7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000323129.fasta
153a12bd34ec3c00f6e9504f91c455de94bccea581026841c99f93d3dc538b18  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000401081.fasta
205924ab6e6b5b88ccaa8d913225ef0f0ae3adf2e05e447e98b542b3fd61f2e7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000504634.fasta
45e2feeca5d2b43235c02e0eb1160737e2174e5319bcdff6936e523b0364bc89  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000649328.fasta
ba24e1285ed71ea1ffea57e4d4313ecea44e888ee84287fc284d3404c71234d1  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000681084.fasta
b1d0c4dcba234929cda8609f97c1698cec05a4239d9de9ce1b13d3b3c6487f62  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000774326.fasta
69eb75d7f49d3a48511f666ff39b0eedbd36b240128c7f06974854a31211c7cf  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000802098.fasta
4a8dc87bf9f6ffb5a1586f5265d9f3cda75ab3ff26ae8e90b17e4e48d21cb615  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP000830882.fasta
aa12133697d12f6a0c81b85f3123ce2c01cf8da875a062b1bd28eefed29c4c82  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001140206.fasta
507df1eb395e7508c0e748f731bc68d4d293b904739fa5696984fc74863da09c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001168540.fasta
5bfcde8b4fe94add6c442c65026b080284f52b47fed0d82d4c2c1e075a6cf1e2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001201262.fasta
a0cd9ee4beb66ecdc4597ab4ddeb01bfe9f37069d59a2f46b7cd7393716e61ff  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001377567.fasta
eb3da5fc595960698f6bfe64aa6afc1d53f19342c9240f579fb034e7fe7ec0a3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001378960.fasta
e05ac10236229dcd2d3ac40518ea562c2b36a57f0b11432c942d89438e5d63e2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001515641.fasta
1e2796849a6c0ae1bc6f46c4689f29e797f5686e37d52f8b181029ddd7e39533  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001519887.fasta
cdaef37ebf5682ceb1464adc626de12c054bc69161778de1f878e77a5265c574  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001610563.fasta
678b50eb97bf144647a8720391ffd1df9f6d8680fb3706338fa9ee2728f8a8a1  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/UP001732720.fasta
136ff88e0c05494f6c9cd09a8c1e0129076b6e0802276601a7940c6e9a5e1d56  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Acanthochromis_polyacanthus.fasta
95b6bc4a0134a9054351548ef2e0935ee6fd45b54b5f9ec41856c6802923eba7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Aeromonas_veronii.fasta
4f33e78803ced5114584d5f1c7a0ce8cdf83c5534596eef5e61a4eb3f7ffc348  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Arabidopsis_thaliana.fasta
90be1ce93705f65517cfe3fbc9789648debe7abf22ee4bc11bf7aad27fdd1ed0  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Bos_taurus.fasta
ba24e1285ed71ea1ffea57e4d4313ecea44e888ee84287fc284d3404c71234d1  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Bothrops_jararaca.fasta
56c85cdd4cca5663f7e439a666769016aaec05bf733f7d189ae67b6218faff5c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Botrytis_cinerea.fasta
4e1ca0066dd2fb34c12d10f62b4cdd4ff3bce6415994bb9fb053b73c75af0793  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Caenorhabditis_elegans.fasta
5e49f15913722d518ddc7662fb2fd663de6ac72e8200c1dd7c202b70e79694e9  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Candida_albicans.fasta
e727da43fa35a3f384537cf9324aaaf27654aa3101ce4293fd53839d3c48f387  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Chlamydomonas_reinhardtii.fasta
8781104df9e8dd1cdf53b89292e855c9af272d48291dfba10747ab6b0cee309c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Ciona_intestinalis.fasta
3d63b778aa8bfeae49d5a8bc5096ce18707229e1a9049341a720ed9d27f2712a  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Danio_rerio.fasta
751bffc6f9718e75f41b25396b937814229169efa87cd5f47018e9075698708d  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Drosophila_melanogaster.fasta
20fc109d91ed4b834bd26cc1fc04d727b6fb04ec2ed28ef2ffd81b9dc15c0e47  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Escherichia_coli.fasta
dbf259de0c637fb24010152bdf7249af733d3e4830ded9e953d39bc679edd0df  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Gallus_gallus.fasta
83aa9913c0c86de19d7cb3f786bb918fe3249ebc3d2a8a2868e5c01f6c2694ca  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Geobacillus_kaustophilus.fasta
6f649d9224beabad4285cbed762455367f5d0387c0862bacd7314449c517f417  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Homo_sapiens.fasta
795d1caf535c43d007fdac3630ac5c5897016d9d7e73212913c20522493c11f2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Hordeum_vulgare.fasta
d98bc3865d6f8a76314ec098774e7deea323484aa36f92390f5191ff86eac88d  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Listeria_monocytogenes.fasta
195c62b3039fae76746d1d0556b5f2e9bbc03294511ab3010450b857013dfcea  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Macaca_fascicularis.fasta
d187aa22b8cd239d730867746b038f9ed154f6063a4a552a5472d90cafa48198  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Marchantia_polymorpha.fasta
4f3e21581a1b00173ff992f9ebba26a6918f964af1ca67d073551acc96c87c9c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Mus_musculus.fasta
ca378270619c48caa9722d83b30f305bf1d7785df57d7a40852a97a036e307ac  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Neisseria_gonorrhoeae.fasta
c640d77b8237453d5df0c2a94eccdd10eca76cae7629e2f40c40bd9c250c1519  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Oryctolagus_cuniculus.fasta
eb3da5fc595960698f6bfe64aa6afc1d53f19342c9240f579fb034e7fe7ec0a3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Pichia_kluyveri.fasta
ee43e442f799fdfef875ef7c4c8a6b7c10edc9e4b63208c590fa63389f94dac3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Pichia_kudriavzevii.fasta
03cfbfaa623e4969186bea40d2f10e0a86ca44a1c06b90a5717d5d436838ecc2  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Plasmodium_falciparum.fasta
65108755b1afec131546a09cfe8ebb6168a143c9e9f4c04497620a6c5bce4b5e  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Pseudomonas_aeruginosa.fasta
7bfaa4e8d832e2bc09bf5eebac9765b7d386e1e9a3ffd846c47448d42b9c5257  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Rattus_norvegicus.fasta
69eb75d7f49d3a48511f666ff39b0eedbd36b240128c7f06974854a31211c7cf  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Rubrivivax_benzoatilyticus.fasta
ff67ac6e6c9192838c1705199ccd9e9c751f372e22264846e9a39e8756e7a105  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Saccharomyces_cerevisiae.fasta
af263c7cd8ebc5bda1f121ff9b08d748486d0810ee7b1ea2672a50f4c9c216b7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Salmonella_enterica.fasta
af263c7cd8ebc5bda1f121ff9b08d748486d0810ee7b1ea2672a50f4c9c216b7  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Salmonella_typhimurium.fasta
16882ed6dc4b9d472fd0d034be8790556463a149a241b473701a07076c90d436  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Schizosaccharomyces_pombe.fasta
8cf52ef95421dca415e03d3d5a6b698a464b0daf8d1caa1368706fa8be43aef3  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Solanum_lycopersicum.fasta
bf2de2111d2ddd741442afa944171ce65a4d5f4aa84c0913dc8376df6781525a  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Solanum_tuberosum.fasta
55f8c677af2d5ffceeae7f7fabbbf0f94d35c417e6bae1467fdb457a7be9f854  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Staphylococcus_aureus.fasta
bc7b19b38120bcbfdd6ca9548b3f7545dd9602f66f3f0a273d238f9e6f99da62  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Streptococcus_pneumoniae.fasta
0747e88316397f640632153c5a740fbb486074cec72ca4c1b1177f79f664167c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Streptomyces_albidoflavus.fasta
091c88de68f45f65056b9fb128245363f70915bc16e539a08f8ac0adf58c80f8  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Sus_scrofa.fasta
1691716754e1226e946200c67e9c21d9514c167b0a74bd133422cda36d4f6f03  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Triticum_aestivum.fasta
b9cc75c92fe801f89fef36e5726b64b9d771441b53f69504a0c1f64b5c45381c  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Trypanosoma_brucei.fasta
ba52f264befb6eeaf9741078ed2354dddd761751c61e2e24c68da8c411be1e24  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Xenopus_laevis.fasta
5083aeb71fc4851f8a17e3e3903d1ac969ee47ebb87c7cb9db84ce7bbdc829bf  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/panel/Zea_mays.fasta
6d19e9a77e4228348a3c26fe8f764f90c74c66813248659b65bc4718f91dff42  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/db/requests.tsv
4b0e6e97ab1d618baa38be612a787d884a781e07f627b919c4b29ac064db5382  /ceph/ibmi/abi/oliver/prose-auto-20261029/search/db/crap.fasta
```

## Step 2: executables, bundles, environment, caches (2026-11-22 14:14 CET)

```
tree 65b38d42c3c3ee2ab5d3f9bb3e1a06c80f6b810a99c1d1a0747a52bf4cc18812  /ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag  (FASTag 1.15.3)
tree c89419390539de9ba31acaefaf823559a9922806ce80e116fd2d581422e03c35  /ceph/ibmi/abi/oliver/opt/fastag-1.10.0/FASTag  (FASTag 1.10.0 (fast_0110))
tree 2f2c380ceae2a4d540ffe492a7e48ffe2e41b13136239c3abae855441fd4c5f5  /ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr  (ProSE (corpus build))
tree 1d5b917b3eda24806fdbc0bdbe262d87b58ecea628d0af46b633aa79673e6e91  /ceph/ibmi/abi/oliver/opt/prose-check-20261028-deps/lib  (ProSE libraries)
tree ef5838fe0ef01c7fa903aa81bd360c1c6e8c32aeb08d973d62e0791154438ce2  /ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169  (FileConverter (diaweaver build, ddaPASEF))
tree ff89f24aa0d68aab9296ad691ecf427154bb64f23a6d7e6d4ba1b65874da1bd2  /ceph/ibmi/abi/oliver/opt/openms-develop-70cebcd144-vendor  (FileConverter (mirror recipe, Thermo))
tree eecab213bcf91fe0ffe9aa885f0050f3e312ec5b92286a73124ca8497dfe7bff  /ceph/ibmi/abi/oliver/opt/mzpeak-convert/0.18.0  (mzpeak-convert 0.18.0)
tree c4438fed49d4c32755a095a2cd1e3e65ad3d5b249af4c5c572b99131043f2daf  /ceph/ibmi/abi/oliver/prose-auto-20261029/opt/FASTag-tax-ext/share-FASTag-taxonomy  (FASTag taxonomy)
tree c0a726833c79cbe1cb4e2fedd8ce849c9daee3f3a60e68b1f7046de4f5c601cd  /ceph/ibmi/abi/oliver/prose-auto-20261029/prose-auto/prose_auto  (prose_auto 0.1.12)
tree 0ae10988e8ff6e408394d9089f57061a02a047889c59ae30ea940ffa730efb99  /ceph/ibmi/abi/oliver/prose-auto-20261029/harness/v0110c_code  (prose_auto 0.1.10 (fast_0110))
tree 6143e5a4e6f769fe5c10a0b5b9a3493bfb573e9743669e371445b3e227cce8f0  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/jev_cache  (Jev cache)
tree 82040934fccd588f4f768f975d145a484365cbf2c3ef99cace7ff0e29d786b34  /ceph/ibmi/abi/oliver/prose-auto-20261029/text_cache_v019  (text cache)
bb96332bf36b03fb2c62a50e26e035cae5965ba9c9705191dbf598e5b7b3b5a6  /ceph/ibmi/abi/oliver/opt/fastag-1.15.3/FASTag/FASTag
bb96332bf36b03fb2c62a50e26e035cae5965ba9c9705191dbf598e5b7b3b5a6  /ceph/ibmi/abi/oliver/opt/fastag-1.10.0/FASTag/FASTag
60ed9c301485139165500f8643d19d174b934ca244964817041c7c31c981eda7  /ceph/ibmi/abi/oliver/opt/prose-check-20261028/pr/bin/ProSE
beaae92b8f07c501a47c4037239c75ee153752850fe18cf68508ad47511578b6  /ceph/ibmi/abi/oliver/spextract/diaweaver-7ed1169/bin/FileConverter
248f45f323057474767ca3f102bfe166b743e636a070170c27243ef514e5001a  /ceph/ibmi/abi/oliver/opt/openms-develop-70cebcd144-vendor/bin/FileConverter
0b4d1a0d1e29c81c2849b956515ecc91828068836d96e8ac7eff2b1395d50798  /ceph/ibmi/abi/oliver/opt/mzpeak-convert/0.18.0/mzpeak-convert
7781287e9b1de155dbe572df16e014491f01983c3738066405517fd585106dd0  /ceph/ibmi/abi/oliver/prose-auto-20261029/opt/bin/diaspextract
4b0e6e97ab1d618baa38be612a787d884a781e07f627b919c4b29ac064db5382  /ceph/ibmi/abi/oliver/prose-auto-20261029/search/db/crap.fasta
d259efe7fef73e3af4ffe5aee474b1f40cad0b31b9b7155a8cecdaafca080d4e  /ceph/ibmi/abi/oliver/prose-auto-20261029/entrap.py
3b714d55680d34285caaa77671c800928393d6cf917d3f4127392284ab3bd5e2  /ceph/ibmi/abi/oliver/prose-auto-20261029/corpus_ids/corpus.py
40a0daa88845190d6805a8c5f99ebc22cc998feb567e95814d7a059dd3bb9ba6  /home/kohlbach/.claude/skills/pride/tools.json
670a2d1b25648f9d69000b25e3e2cf01adb4a5e985e5e60fbf268366d30feb1a  /ceph/ibmi/abi/oliver/prose-auto-20261029/v0112/heldout/pip_freeze.txt
```

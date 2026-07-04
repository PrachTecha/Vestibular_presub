"""Shared gene sets for the vulnerability-index analysis (Fig 2i-l).

`metabolic_load`:   OXPHOS complexes (I, IV, V) + TCA + ion pumps + Ser/Gly/Gln
                    metabolism + ribosomal + splicing factors / helicases.
                    Represents biosynthetic / translational / energetic demand.

`proteostasis_capacity`: Chaperones, proteasome, ubiquitin, autophagy,
                         lysosomal hydrolases, ERAD / UPR. Represents the
                         capacity to fold, refold, degrade, or recycle proteins.

Vulnerability index = log2(metabolic_load + eps) - log2(proteostasis_capacity + eps)
                      (interpretation: high biosynthetic demand vs low coping
                      capacity = vulnerable cell)
"""

metabolic_genes = [
    # Complex I (Nduf*)
    "Ndufa1", "Ndufa2", "Ndufa3", "Ndufa4", "Ndufa5", "Ndufa6", "Ndufa7",
    "Ndufa8", "Ndufa9", "Ndufa10", "Ndufa11", "Ndufa12", "Ndufa13",
    "Ndufb1", "Ndufb2", "Ndufb3", "Ndufb4", "Ndufb5", "Ndufb6", "Ndufb7",
    "Ndufb8", "Ndufb9", "Ndufb10", "Ndufb11", "Ndufc1", "Ndufc2",
    "Ndufs1", "Ndufs2", "Ndufs3", "Ndufs4", "Ndufs5", "Ndufs6", "Ndufs7", "Ndufs8",
    # Complex IV (Cox*)
    "Cox4i1", "Cox5a", "Cox5b", "Cox6a1", "Cox6b1", "Cox6c",
    "Cox7a2l", "Cox7b", "Cox7c",
    # Complex V (Atp5*)
    "Atp5f1", "Atp5g1", "Atp5g2", "Atp5h", "Atp5i", "Atp5j", "Atp5j2",
    "Atp5k", "Atp5l", "Atp5o",
    # TCA
    "Aco2", "Idh3a", "Idh3b", "Idh3g", "Ogdh", "Dlst",
    "Suclg1", "Sucla2", "Suclg2", "Sdha", "Sdhb", "Sdhc", "Sdhd", "Mdh2", "Fh1",
    # Ion pumps (Na/K-ATPase, Ca-ATPase)
    "Atp1a1", "Atp1a2", "Atp1b1", "Atp1b2", "Atp2a2", "Atp2b1", "Atp2b2",
    # Amino-acid metabolism (Ser/Gly/Gln)
    "Psat1", "Psph", "Phgdh", "Shmt2", "Asns", "Gls", "Gpt2", "Got2",
    # Ribosomal large subunit (Rpl*)
    "Rpl3", "Rpl5", "Rpl7", "Rpl10", "Rpl11", "Rpl23", "Rpl26", "Rpl27",
    "Rpl28", "Rpl29", "Rpl30", "Rpl31", "Rpl32", "Rpl34", "Rpl35", "Rpl36",
    "Rpl37", "Rpl38",
    # Ribosomal small subunit (Rps*)
    "Rps3", "Rps5", "Rps8", "Rps10", "Rps11", "Rps12", "Rps13", "Rps14",
    "Rps15", "Rps16", "Rps17", "Rps18", "Rps19", "Rps20", "Rps21", "Rps23",
    "Rps24", "Rps25", "Rps26", "Rps27", "Rps28",
    # Splicing factors
    "Hnrnpa1", "Hnrnpa2b1", "Hnrnpd", "Hnrnph1", "Hnrnpm",
    "Srsf1", "Srsf2", "Srsf3", "Srsf7",
    "Snrpa1", "Snrpb", "Snrpd1", "Snrpd2", "Snrpe", "Snrpf", "Snrpg",
    # DEAD-box helicases
    "Ddx5", "Ddx17", "Ddx39b",
]

proteostasis_genes = [
    # Chaperones (HSP70 / HSP90 / small HSP)
    "Hspa1a", "Hspa1b", "Hspa4", "Hspa5", "Hspa8",
    "Hsp90aa1", "Hsp90ab1", "Hspb1",
    # TRiC / CCT chaperonin
    "Cct1", "Cct2", "Cct3", "Cct4", "Cct5", "Cct6a", "Cct7", "Cct8", "Tcp1",
    # 20S proteasome (alpha / beta)
    "Psma1", "Psma2", "Psma3", "Psma4", "Psma5", "Psma6", "Psma7",
    "Psmb1", "Psmb2", "Psmb3", "Psmb4", "Psmb5", "Psmb6", "Psmb7", "Psmb8",
    # 19S regulatory (ATPase + non-ATPase)
    "Psmc1", "Psmc2", "Psmc3", "Psmc4", "Psmc5", "Psmc6",
    "Psmd1", "Psmd2", "Psmd3", "Psmd4", "Psmd5", "Psmd6", "Psmd7",
    "Psmd8", "Psmd9", "Psmd10", "Psmd11", "Psmd12", "Psmd13", "Psmd14",
    # Ubiquitin-conjugating
    "Uba1", "Ube2d1", "Ube2d2", "Ube2d3", "Ube2n", "Ube2k",
    # Autophagy core
    "Atg3", "Atg5", "Atg7", "Atg10", "Atg12", "Atg16l1",
    "Becn1", "Map1lc3a", "Map1lc3b",
    # Lysosomal hydrolases / membrane
    "Ctsb", "Ctsd", "Ctsh", "Ctsk", "Ctss", "Lamp1", "Lamp2",
    # ERAD / UPR
    "Syvn1", "Derl1", "Derl2", "Herpud1",
    "Ern1", "Atf6", "Xbp1", "Ddit3",
]

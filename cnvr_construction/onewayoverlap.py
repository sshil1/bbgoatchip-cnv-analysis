import pandas as pd

# Load your validated file
df = pd.read_csv("/home/gmbluser/803m_concordance_validated.tsv", sep='\t')

# Filter for short-read calls only
sr = df[df['source'] == 'cnvpytor'].copy()

# Since your previous script already marked 'matched' = 'yes' based on overlap
# Let's look at the DELs specifically as they are your most reliable set
dels = sr[sr['svtype'] == 'DEL']
dups = sr[sr['svtype'] == 'DUP']
mixed = sr[sr['svtype'] == 'MIXED']

print("--- RELIABILITY TIER LIST ---")
print(f"Deletions:    {len(dels[dels['matched']=='yes'])/len(dels)*100:.1f}% Reliable")
print(f"Mixed:        {len(mixed[mixed['matched']=='yes'])/len(mixed)*100:.1f}% Reliable")
print(f"Duplications: {len(dups[dups['matched']=='yes'])/len(dups)*100:.1f}% Reliable")

# ACTIONABLE: Save a 'Golden Set' of validated variants
golden_set = sr[sr['matched'] == 'yes']
golden_set.to_csv("803m_validated_high_confidence.tsv", sep='\t', index=False)
print(f"\nSaved {len(golden_set)} high-confidence variants to 803m_validated_high_confidence.tsv")

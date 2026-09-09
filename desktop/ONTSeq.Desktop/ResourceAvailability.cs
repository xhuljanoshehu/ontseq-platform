namespace ONTSeq.Desktop;

public enum ResourceFamilyAvailability
{
    Ready,
    NotInstalled,
    Incomplete,
    Unavailable
}

public sealed record ResourceFamilyState(
    string GenomeBuild,
    ResourceFamilyAvailability Availability,
    string Detail)
{
    public bool CanAnalyze => Availability == ResourceFamilyAvailability.Ready;
    public bool CanInstall => Availability == ResourceFamilyAvailability.NotInstalled;
    public bool CanRepair =>
        Availability is ResourceFamilyAvailability.Ready or ResourceFamilyAvailability.Incomplete;

    public string AvailabilityLabel => Availability switch
    {
        ResourceFamilyAvailability.Ready => "lokal bereit",
        ResourceFamilyAvailability.NotInstalled => "nicht installiert",
        ResourceFamilyAvailability.Incomplete => "unvollständig",
        _ => "Status nicht verfügbar"
    };
}

public sealed record DesktopProfileOption(
    DesktopAnalysisProfile Profile,
    ResourceFamilyState ResourceFamily)
{
    public bool IsEnabled => ResourceFamily.CanAnalyze;
    public string AvailabilityLabel => ResourceFamily.AvailabilityLabel;

    public override string ToString() => $"{Profile.DisplayName} · {AvailabilityLabel}";
}

public sealed record ProfileSelectionDecision(
    DesktopAnalysisProfile? Profile,
    bool IsFallback,
    string? Notice);

public static class DesktopResourcePolicy
{
    public static IReadOnlyList<string> GenomeBuilds { get; } =
        Array.AsReadOnly(new[] { "GRCh37", "GRCh38" });

    public static IReadOnlyDictionary<string, ResourceFamilyState> UnavailableFamilies(
        string detail) => GenomeBuilds.ToDictionary(
        build => build,
        build => new ResourceFamilyState(
            build,
            ResourceFamilyAvailability.Unavailable,
            detail),
        StringComparer.Ordinal);

    public static IReadOnlyList<DesktopProfileOption> BuildProfileOptions(
        IReadOnlyDictionary<string, ResourceFamilyState> families) =>
        DesktopProfiles.Supported.Select(profile =>
        {
            var state = families.TryGetValue(profile.GenomeBuild, out var available)
                ? available
                : new ResourceFamilyState(
                    profile.GenomeBuild,
                    ResourceFamilyAvailability.Unavailable,
                    "Für diesen Build wurde kein Ressourcenstatus geliefert.");
            return new DesktopProfileOption(profile, state);
        }).ToArray();

    public static ProfileSelectionDecision ResolveInitialProfile(
        string? configuredProfileId,
        IReadOnlyDictionary<string, ResourceFamilyState> families)
    {
        var configured = DesktopProfiles.Supported.FirstOrDefault(profile => string.Equals(
            profile.ProfileId, configuredProfileId, StringComparison.Ordinal));
        if (configured is not null && IsReady(configured, families))
            return new ProfileSelectionDecision(configured, false, null);

        var fallback = DesktopProfiles.Supported.FirstOrDefault(profile => IsReady(profile, families));
        if (fallback is null)
        {
            return new ProfileSelectionDecision(
                null,
                false,
                "Im aktuellen Resource-Root ist noch keine Build-Familie lokal bereit.");
        }

        var configuredLabel = configured?.ProfileId ?? configuredProfileId ?? "kein Profil";
        return new ProfileSelectionDecision(
            fallback,
            true,
            $"{configuredLabel} ist im aktuellen Resource-Root nicht bereit; " +
            $"{fallback.ProfileId} wurde als erstes lokal bereites Profil gewählt.");
    }

    public static bool CanPersistProfile(
        string profileId,
        IReadOnlyDictionary<string, ResourceFamilyState> families)
    {
        var profile = DesktopProfiles.Supported.FirstOrDefault(item => string.Equals(
            item.ProfileId, profileId, StringComparison.Ordinal));
        return profile is not null && IsReady(profile, families);
    }

    private static bool IsReady(
        DesktopAnalysisProfile profile,
        IReadOnlyDictionary<string, ResourceFamilyState> families) =>
        families.TryGetValue(profile.GenomeBuild, out var state) && state.CanAnalyze;
}

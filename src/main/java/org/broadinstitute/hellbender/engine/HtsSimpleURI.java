package org.broadinstitute.hellbender.engine;

import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.Utils;

import java.net.URI;
import java.net.URISyntaxException;
import java.nio.file.*;
import java.nio.file.spi.FileSystemProvider;

/**
 * Default implementation for HtsURI.
 *
 * This class takes raw strings that are to be interpreted as URIs, and converts them internally to a URI or
 * Path object. If no scheme is provided as part of the raw string used in the constructor(s), the URI is assumed
 * to represent a file on the local file system, and will be backed by a URI with a "file:/" scheme and a path
 * part that is automatically encoded/escaped to ensure it is a valid URI. If the raw string contains a scheme,
 * it will be backed by a URI formed from the raw string as presented, with no further encoding/escaping.
 *
 * For example, a URI that contains a scheme, and has an embedded "#" in the path will be treated as having a
 * fragment delimiter. If the URI contains an scheme, the "#" will be escaped.
 *
 * There are 3 succeeding levels of validation:
 *
 * 1) HtsURI constructor - requires a syntactically valid URI, possibly containing a scheme (if no scheme is present
 * the path part will be escaped/encoded)
 * 2) isNio - returns true if the URI is syntactically valid, and there is an installed NIO provider that matches
 * the URI scheme
 * 3) isPath - syntactically valid URI that can be resolved to a java.io.Path by the provider
 *
 * <scheme>:<scheme-specific-part>
 * <scheme>://<authority><path>?<query>
 *  absoluteURI   = scheme ":" ( hier_part | opaque_part )
 *       hier_part     = ( net_path | abs_path ) [ "?" query ]
 *       net_path      = "//" authority [ abs_path ]
 *       abs_path      = "/"  path_segments
 *
 * A URI is absolute if, and only if, it has a scheme component.
 * A URI is opaque if, and only if, it is absolute (has a scheme) and its
 * scheme-specific part does not begin with a slash character ('/')
 *
 * A relative reference that does not begin with a scheme name or a
 * slash character is termed a relative-path reference.
 *
 * A relative reference beginning with a single slash character is
 * termed an absolute-path reference, as defined by <abs_path> in
 * Section 3.
 *
 * URI that do not make use of the slash "/" character for separating
 * hierarchical components are considered opaque by the generic URI
 * parser.
 *
 * TODO: for simplicity, this currently contains some (cloud) functionality and tests that should ultimately
 * TODO: be moved into GatkURI, but its easier for now to see everything together
 * TODO: should we use the normalized uri string for hash/compare purposes, or the user-supplied string ?
 * TODO: special handling for authority/credentials, query params, fragment?
 *
 * TODO: fix FeatureInput.makeIntoAbsolutePath - makes a File Path out of a gendb path. amongst other things
 */
public class HtsSimpleURI implements HtsURI {

    private final String    rawInputString;     // raw input string provided by th user; may or may not have a scheme
    private final URI       uri;                // working URI; always has a scheme (assume "file" if not provided)
    private Path            cachedPath;         // cache the Path associated with this URI if its "Path-able"
    private String          pathFailureReason;  // cache the reason for "toPath" conversion failure

    public HtsSimpleURI(final String rawInputString) {
        Utils.nonNull(rawInputString);
        this.rawInputString = rawInputString;

        URI tempURI;
        try {
            // If the input URI already has a scheme (including a "file" scheme), we assume its already properly
            // escaped. If no scheme component is present, then assume its a raw path on the local file system, so
            // try to get a Path first, and then recreate the URI by retrieving it from the resulting Path. This
            // ensures that input strings that contain embedded characters that would otherwise be interpreted as
            // URI syntax (such as embedded fragment specifiers ("#") that are valid in file names) are properly
            // escaped.
            tempURI = new URI(rawInputString);
            if (!tempURI.isAbsolute()) {
                // NOTE: This case (no scheme) is the only one where we resolve the URI to a Path at construction time.
                try {
                    setCachedPath(Paths.get(rawInputString));
                    tempURI = getCachedPath().toUri();
                } catch (InvalidPathException e) {
                    throw new IllegalArgumentException(e.getMessage(), e);
                }
            }
            uri = tempURI;
            if (!uri.isAbsolute()) {
                // assert the invariant that every URI we create has a scheme, even if the raw input string does not
                throw new GATKException("URI has no scheme");
            }
            //displayComponents(uri);
        } catch (URISyntaxException e) {
            final String errorMessage = String.format("%s must be a valid URI. '%s'/'%s'", rawInputString, e.getMessage(), e.getReason());
            throw new IllegalArgumentException(errorMessage);
        }
    }

    /**
     * Converts the URI to a {@link Path} object.
     *
     * @return the resulting {@code Path}
     * @throws UserException if an I/O error occurs when creating the file system
     */
    @Override
    public Path toPath() {
        if (getCachedPath() != null) {
            return getCachedPath();
        } else {
            final Path tmpPath = Paths.get(getURI());
            setCachedPath(tmpPath);
            return tmpPath;
        }
    }

    @Override
    public URI getURI() {
        return uri;
    }

    @Override
    public String getURIString() {
        return rawInputString;
    }

    @Override
    public boolean isPath() {
        try {
            return getCachedPath() != null || toPath() != null;
        } catch (ProviderNotFoundException |
                FileSystemNotFoundException |
                IllegalArgumentException |
                UserException |
                AssertionError e) {
            // jimfs throws an AssertionError that wraps a URISyntaxException when trying to create path where
            // the scheme-specific part is missing or incorrect
            pathFailureReason = e.getMessage();
            return false;
        }
    }

    @Override
    public boolean isNIO() {
        // try to find a provider; assume that our URI always has a scheme
        for (FileSystemProvider provider: FileSystemProvider.installedProviders()) {
            if (provider.getScheme().equalsIgnoreCase(uri.getScheme())) {
                return true;
            }
        }
        return false;
    }

    @Override
    public String getToPathFailureReason() {
        if (pathFailureReason == null) {
            try {
                toPath();
                return String.format("'%s' appears to be a valid Path", rawInputString);
            } catch (ProviderNotFoundException e) {
                return String.format("ProviderNotFoundException: %s", e.getMessage());
            } catch (FileSystemNotFoundException e) {
                return String.format("FileSystemNotFoundException: %s", e.getMessage());
            } catch (IllegalArgumentException e) {
                return String.format("IllegalArgumentException: %s", e.getMessage());
            } catch (UserException e) {
                return String.format("UserException: %s", e.getMessage());
            }
        }
        return pathFailureReason;
    }

    // get the cached path associated with this URI if its already been created
    protected Path getCachedPath() { return cachedPath; }

    protected void setCachedPath(Path path) {
        this.cachedPath = path;
    }

    @Override
    public String toString() {
        return rawInputString;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof HtsSimpleURI)) return false;

        HtsSimpleURI that = (HtsSimpleURI) o;

        if (!getURIString().equals(that.getURIString())) return false;
        if (!getURI().equals(that.getURI())) return false;
        return true;
    }

    @Override
    public int hashCode() {
        int result = getURIString().hashCode();
        result = 31 * result + getURI().hashCode();
        return result;
    }

    private void displayComponents(final URI uri) {
        System.out.println(
                String.format(
                        "Original:         %s\nURI:              %s\nASCII:            %s\nNormalized:       %s\nNormalized ASCII: %s\n",
                        rawInputString,
                        uri.toString(),
                        uri.toASCIIString(),
                        uri.normalize().toString(),
                        uri.normalize().toASCIIString()
                )
        );

        System.out.println(
                String.format(
                        "Scheme:         %s\nAuthority:      %s\nHost:           %s\n\n",
                        uri.getScheme(),
                        uri.getAuthority(),
                        uri.getHost()
                )
        );
    }

}
